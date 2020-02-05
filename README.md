SlowerLogLog
==

A sample implementation of the cardinality estimation algorithm described here:

https://www.evanmiller.org/slower-log-log.html

This repository contains a Julia implementation (SlowerLogLog.jl) as well as a
C command-line tool, which can be piped input.

Compilation:

```
clang slowcount.c -o slowcount
```

Example usage:

```
ls /tmp/ | ./slowcount
```
