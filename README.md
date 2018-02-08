# Stripeline2

Next generation of the STRIP pipeline, written in [Julia](https://julialang.org) (see [here](https://github.com/lspestrip/stipeline) for the previous version, written in Python).

Early in 2018, we decided to rewrite the STRIP simulation pipeline in Julia because of the following issues with the old version:

1. Significant amount of memory was required
2. Code was quite verbose
3. In order to improve the execution speed, we had to implement a few routines in Fortran

This new version of the pipeline should require less memory
and perform somewhat faster, although full-scale benchmarks
have yet to be done.
