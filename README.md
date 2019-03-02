# Stripeline.jl

[![Build Status](https://travis-ci.org/lspestrip/Stripeline.jl.svg?branch=master)](https://travis-ci.org/lspestrip/Stripeline.jl)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://lspestrip.github.io/Stripeline.jl/latest)
[![Coverage Status](https://img.shields.io/coveralls/lspestrip/Stripeline.jl.svg)](https://coveralls.io/r/lspestrip/Stripeline.jl?branch=master)

Next generation of the STRIP pipeline, written in [Julia](https://julialang.org)
(see [here](https://github.com/lspestrip/stipelinepy) for the previous version,
written in Python).

Early in 2018, we decided to rewrite the STRIP simulation pipeline in Julia
because of the following issues with the old version:

1. Significant amount of memory was required
2. Code was quite verbose
3. In order to improve the execution speed, we had to implement a few
   routines in Fortran: this meant that the code was harder to read,
   as knowledge of two languages (Fortran+Python) was required.

This new version of the pipeline should require less memory and perform somewhat
faster, although full-scale benchmarks have yet to be done. It is far more
advanced than the older pipeline, as the instrument database is updated and a
destriper is already implemented.

## How to contribute

See the file
(https://github.com/lspestrip/Stripeline.jl/blob/devel/CONTRIBUTING.md)[CONTRIBUTING.md].
