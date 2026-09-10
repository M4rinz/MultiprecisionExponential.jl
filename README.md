# MultiprecisionExponential

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://M4rinz.github.io/MultiprecisionExponential.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://M4rinz.github.io/MultiprecisionExponential.jl/dev/)
[![Build Status](https://github.com/M4rinz/MultiprecisionExponential.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/M4rinz/MultiprecisionExponential.jl/actions/workflows/CI.yml?query=branch%3Amain)

**Arbitrary precision scaling and squaring algorithm for the exponential of a matrix**

`MultiprecisionExponential.jl` implements the algorithm by Fasi and Higham to compute the exponential of a matrix at arbitrary precision.

# Compatibility
Julia v1.12 or higher (maybe lower as well).

# Example
```julia
using MultiprecisionExponential

A = rand(2,2)

Y = exp_mp(A, working_precision=256)
```

# How to install
Run Julia, then... Wait for me to register the package...

# References

**Article:**
> M. Fasi and N. J. Higham, An Arbitrary Precision Scaling and Squaring Algorithm for the Matrix Exponential
> SIAM J. Matrix Anal. Appl., Vol. 40.4 (2019), pp.1233-1256.
> [doi: 10.1137/18M1228876](https://doi.org/10.1137/18M1228876)

**Original MATLAB implementation** at [this link](https://github.com/mfasi/mpexpm).