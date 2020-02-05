/* C command-line tool for cardinality estimation */
/* Uses an MLE-based variant of HyperLogLog */

#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <math.h>

/* Borrowed from djb2 */
uint32_t
hash_string(unsigned char *str)
{
    uint32_t hash = 5381;
    int c;

    while ((c = *str++))
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */

    return hash;
}

/* Borrowed from Julia */
static inline uint32_t hash_int(uint32_t a) {
    a = a + 0x7ed55d16 + (a << 12);
    a = a ^ 0xc761c23c ^ (a >> 19);
    a = a + 0x165667b1 + (a << 5);
    a = a + 0xd3a2646c ^ (a << 9);
    a = a + 0xfd7046c5 + (a << 3);
    a = a ^ 0xb55a4f09 ^ (a >> 16);
    return a;
}

void update_registers(unsigned char *registers, size_t len, uint32_t value) {
    int i;
    for (i=0; i<len; i++) {
        int k = __builtin_ctz(~value);
        if (k > registers[i])
            registers[i] = k;
        value = hash_int(value);
    }
}

/* Probability of n given register value k */
double g0f(double n, int k) {
    double a = log1p(scalbn(-1.0, -(k+1)));
    double b = log1p(scalbn(-1.0, -k));
    return exp(n * a) - exp(n * b);
}

/* First derivative */
double g1f(double n, int k) {
    double a = log1p(scalbn(-1.0, -(k+1)));
    double b = log1p(scalbn(-1.0, -k));
    return exp(n * a) * a - exp(n * b) * b;
}

/* Second derivative */
double g2f(double n, int k) {
    double a = log1p(scalbn(-1.0, -(k+1)));
    double b = log1p(scalbn(-1.0, -k));
    return exp(n * a) * a * a - exp(n * b) * b * b;
}

/* Initial guess for Newton's method */
double harmonic_mean(unsigned char *registers, size_t len) {
    double total = 0.0;
    int i;
    for (i=0; i<len; i++) {
        total += scalbn(1.0, -registers[i]);
    }
    return len / total;
}

/* First derivative of the log-likelihood function */
double jacobian(unsigned char *registers, size_t len, double n) {
    double total = 0.0;
    int i;
    for (i=0; i<len; i++) {
        if (registers[i])
            total += g1f(n, registers[i]) / g0f(n, registers[i]);
    }
    return total;
}

/* Second deriviative of the log-likelihood function */
double hessian(unsigned char *registers, size_t len, double n) {
    double total = 0.0;
    int i;
    for (i=0; i<len; i++) {
        double g0 = g0f(n, registers[i]);
        double g1 = g1f(n, registers[i]);
        double g2 = g2f(n, registers[i]);
        if (registers[i])
            total += ( g0 * g2 - g1 * g1 ) / ( g0 * g0 );
    }
    return total;
}

/* The main event */
double estimate(unsigned char *registers, size_t len) {
    double n = harmonic_mean(registers, len);
    int i;
    for (i=0; i<8; i++) {
        n -= jacobian(registers, len, n) / hessian(registers, len, n);
    }
    return n;
}

void usage() {
    printf("Usage: slowcount [ <# registers, between 10 and 16000> ]\n");
}

int main(int argc, char **argv) {
    size_t len = 2000;
    if (argc == 2) {
        len = strtol(argv[1], NULL, 10);
        if (len < 10 || len > 16000) {
            usage();
            return 1;
        }
    } else if (argc != 1) {
        usage();
        return 1;
    }
    unsigned char *registers  = calloc(len, 1);
    char *line = NULL;
    size_t linecap = 0;
    ssize_t linelen;
    int count = 0;
    while ((linelen = getline(&line, &linecap, stdin)) > 0) {
        uint32_t hash = hash_int(hash_string((unsigned char *)line));
        update_registers(registers, len, hash);
    }
    double est = estimate(registers, len);
    double var = -1.0/hessian(registers, len, est);

    free(registers);

    printf("%.2lf ± %.2lf\n", est, sqrt(var));
    return 0;
}
