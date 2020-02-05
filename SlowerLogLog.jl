#!/usr/bin/julia

using Printf

function g(d, n, k)
    if k == 0
        return 0
    end

    return (1 - 0.5^(k+1))^n * (log(1-0.5^(k+1)))^d - (1 - 0.5^k)^n * (log(1-0.5^k))^d
end

function g0(n, k)
    return g(0, n, k)
end

function g1(n, k)
    return g(1, n, k)
end

function g2(n, k)
    return g(2, n, k)
end

function jacobian(registers, n)
    total = 0.0
    for k in registers
        if k > 0
            total += g1(n, k) / g0(n, k)
        end
    end
    return total
end

function hessian(registers, n)
    total = 0.0
    for k in registers
        if k > 0
            total += ( g0(n, k) * g2(n, k) - (g1(n, k))^2 ) / (g0(n, k)^2)
        end
    end
    return total
end

function harmonic_mean(registers)
    total = 0.0
    for k in registers
        total += 0.5^k
    end
    return 1/total * length(registers)
end

function estimate(registers)
    n = harmonic_mean(registers)

    for j = 1:8
        n -= jacobian(registers, n) / hessian(registers, n)
    end

    return n
end

registers = zeros(2000)

# values = 1000:10000

for v in readlines()
    h = hash(v)
    for i = 1:length(registers)
        k = trailing_zeros(h)
        if k > registers[i]
            registers[i] = k
        end
        h = hash(h)
    end
end

n = estimate(registers)

variance = -1/hessian(registers, n)

@printf("%.2lf ± %.2lf\n", n, sqrt(variance))
#println("Harmonic Mean Estimate: ", harmonic_mean(registers))
#println("MLE Estimate: ", n)
#println("Std. error: ", sqrt(variance))
#println("RSE: ", sqrt(variance)/n)
#println("RSE * sqrt(M): ", sqrt(variance)/n*sqrt(length(registers)))
