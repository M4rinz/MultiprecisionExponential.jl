```@meta
CurrentModule = MultiprecisionExponential
```

# MultiprecisionExponential

Documentation for [MultiprecisionExponential](https://github.com/M4rinz/MultiprecisionExponential.jl).


## Installation
To install `MultiprecisionExponential.jl`... well... Wait...

## Example
```julia
using MultiprecisionExponential

A = rand(2,2)

Y_53 = exp_mp(A)
Y256 = exp_mp(A, working_precision=256)
```

## References
Original article, and original MATLAB implementation
