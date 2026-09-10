module MultiprecisionExponential

using LinearAlgebra
using GenericSchur
using LinearMaps, MatrixEquations

const VALID_APPRX = (:taylor, :diagonalcheap)
const VALID_ALGS  = (:transfree, :realschur, :complexschur)


include("structs.jl")
include("subfunctions.jl")
include("exp_mp.jl")

export exp_mp

end
