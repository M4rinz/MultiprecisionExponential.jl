const p_factor = 1.2    # internal increase in precision for the BigFloats


############ Evaluation of matrix polynomials ############
"""
    polyvalm_mp!(Apows, s, β; outputclass=nothing)

Evaluates the matrix polynomial ``\\sum_{k=0}^m\\beta_k(2^{-s}A)^k``
using the Paterson-Stockmeyer scheme [ps]. 
The "degree" of the transformed polynomial is ``r=\\lfloor m/ν\\rfloor``
where ``ν=\\lceil\\sqrt{m}\\rceil``.

# Arguments
- `Apows::AandPowsStruct`: A struct containing the powers of the matrix ``A``, 
    *starting from the identity*
- `s::Real`: The scaling factor ``s``
- `β::AbstractVector`: The vector with the polynomial coefficients: 
    ``\\beta=(\\beta_0,\\ldots,\\beta_{m})``
- `outputclass::Type` (Optional. Defaults to `nothing`): The return type. 
    If set to `nothing`, the output will have elements of the same type as `A`. 

# References 
> [^ps] Paterson, Michael S. and Stockmeyer, Larry J., On the Number of Nonscalar Multiplications Necessary to Evaluate Polynomials
> SIAM J. Comp., Vol. 2.1 (1973), pp. 60-66.
> [doi: 10.1137/0202007](https://doi.org/10.1137/0202007)
"""
function polyvalm_ps!(
        Apows::AbstractVector{<:AbstractMatrix},
        s::Real,
        β_vec::AbstractVector;
        outputclass::Union{Type,Nothing}=nothing     
) 
    T = promote_type(eltype.(Apows)...)
    outputclass = something(outputclass, big(T))

    n = LinearAlgebra.checksquare(Apows[2])

    m = length(β_vec) - 1      # degree of the polynomial
    ν = ceil(typeof(m), √m)    # the "batch size"

    if ν == 0
        Y = outputclass(β_vec[1]) * I(n) 
        return Y
    end

    r = fld(m, ν)              # the "degree" of the P.-S. polynomial
    
    for i in length(Apows)+1:ν+1
        if isodd(i)
            i_half = div(i - 1, 2) + 1
            push!(Apows, Apows[i_half]^2)
        else
            push!(Apows, Apows[2] * Apows[i-1])
        end
    end

    mpowers = copy.(Apows)   
    for i in 1:ν+1
        mpowers[i] ./= 2^(2*big(s)*(i-1))
    end
    # convert creates an alias, if Apows already contains matrices of `outputclass`.
    # Thus, the original accuracy of A would be preserved
    mpowers = convert(Vector{AbstractMatrix{outputclass}}, mpowers) 

    # Evaluate the last B-term, the one of degree m - ν*r = m mod ν
    B = setprecision(floor(Int64, p_factor * precision(BigFloat))) do 
        B = outputclass(β_vec[m+1]) * mpowers[m - ν*r + 1]
        for j in m-1:-1:ν*r
            if j == ν*r
                B += outputclass(β_vec[ν*r + 1]) * I(n)
            else
                B += outputclass(β_vec[j+1]) * mpowers[m - ν*r - (m-j) + 1]
            end
        end
        B;
    end

    # Evaluate B-terms of degree ν-1, and evaluate the main polynomial using Horner
    Y = convert(AbstractMatrix{outputclass}, B) # to retain original precision 

    for kk in r-1:-1:0
        # compute B coeff. in slightly higher precision
        B = setprecision(floor(Int64, p_factor * precision(BigFloat))) do 
            B = outputclass(β_vec[ν*kk + 1]) * I(n)  
            for j in 1:ν-1
                B += outputclass(β_vec[ν*kk + j + 1]) * mpowers[j+1]
            end
            B;
        end
        Y = Y * mpowers[ν+1] + convert(AbstractMatrix{outputclass}, B)
    end

    return Y
end


# Paterson-Stockmeyer on Taylor polynomials.
function polyvalm_tay_exp(
    Apows::AbstractVector{<:AbstractMatrix}, 
    m::Integer, s::Real
)
    n = LinearAlgebra.checksquare(Apows[1])

    ν = ceil(typeof(m), √m)    # the "batch size"
    ν == 0 && throw(DomainError(lazy"Polynomial degree m is $(m). A value greater than 1 is expected"))
    r = floor(typeof(m), m/ν)  # the "degree" of the P.-S. polynomial

    scaling = 2^s

    Y = Apows[ν+1] / (scaling * m)
    for j in ν-1:-1:1
        Y = (Y + Apows[j+1]) / (scaling * (m - ν + j))
    end

    Y = Apows[ν+1] * Y;
    for k in r-2:-1:1
        for j in ν:-1:1
            Y = (Y + Apows[j+1]) / (scaling * (k * ν + j))
        end
        Y = Apows[ν+1] * Y
    end

    for j in ν:-1:1
        Y = (Y + Apows[j+1]) / (scaling * j)
    end
    
    Y = Y + I(n)

    return Y
end



############ Approximants to the exponential ############

# Helper to compute even and odd Padé parts
function _pade_even_odd_parts(
    S::AandPowsStruct,
    m::Integer,
    s::Real
)
    c_num = setprecision(floor(Int64, 2*precision(BigFloat))) do
        mb = big(m)
        ranges = 0:mb
        c_num = (factorial(mb)/factorial(2mb)) ./
                (factorial.(mb .- ranges) .* factorial.(ranges)) .*
                factorial.(2mb .- ranges)
        c_num[1] = big(1.)
        c_num
    end

    Uₑ = polyvalm_ps!(S.powers, s, c_num[1:2:end])

    if m ≥ 1
        Uₒ = polyvalm_ps!(S.powers, s, c_num[2:2:end])
        Uₒ *= (S.A / 2^s)
    else
        Uₒ = zero(S.A)
    end

    return Uₑ, Uₒ
end

function eval_pade!(
    S::AandPowsStruct,
    m::Integer,
    s::Real,
    ::Val{:diagonalcheap}
)
    n = LinearAlgebra.checksquare(S.A)

    Uₑ, Uₒ = _pade_even_odd_parts(S, m, s)
    Qₘ = Uₑ - Uₒ
    Rₘ = Qₘ \ (2*Uₒ) + I(n)

    return Rₘ
end

function eval_pade!(
    S::AandPowsStruct,
    m::Integer,
    s::Real,
    ::Val{:taylor}
)
    # `scalar_error_tayl!` has already added the necessary powers
    ν = ceil(typeof(m), √m)
    length(S.powers) < ν+1 && throw(ArgumentError(lazy"S.powers has length $(length(S.powers)), it should be of length $(ν+1)"))

    return polyvalm_tay_exp(S.powers, m, s)
end


"""
    eval_pade!(S, m, s, approximant)

Approximates ``e^(2^{-s}A)`` using the Padé approximant of degree `m`.
The approximant is either the Taylor polynomial ``T_m`` or the diagonalcheap Padé approximant
``r_m``, depending on the value of `approximant`.
The approximant is evaluated using the Paterson-Stockmeyer (an ad-hoc version for the Taylor approximant).
For `:diagonalcheap`, a closed-form formula is used.

# Arguments
- `S::AandPowsStruct`: A struct with the fields `use_taylor`, `A` and `powers`
    - `A::AbstractMatrix`: The matrix whose scaled version we want the exponential of
    - `powers::AbstractVector`: A vector holding the powers of ``A`` if `S.use_taylor` is true, 
                                or ``A^2`` otherwise. See Remark 1 below. 
    - `use_taylor::Bool`: whether the chosen approximant is the Taylor one or not.
- `m::Integer`: The degree of the Padé approximant. *Note*: an optimal 
                degree must be used. See Remark 2 below.
- `s::Real`: The scaling factor.
- `approximant::Symbol`: The approximant to use (`:taylor` or `:diagonalcheap`)

!!! note "Remark 1"
    For the (diagonal) Padé implementation, at least ``I`` and ``A^2`` must be present. Then
    the missing powers will be inserted in `S.powers`. 
    Instead, for the Taylor version, it is assumed that all the required powers 
    are be inside `S.powers`.  This is because when the 
    error bound is evaluated by [`scalar_error_tayl!`](@ref), 
    the missing powers are inserted in `S.powers`. 
    This function should be used in pair with [`scalar_error_tayl!`](@ref), as this code does,
    one must be wary about testing it independently.

!!! note "Remark 2"
    The degree `m` of the Taylor approximant must be one of the "optimal ones" (see article [hf19_mpexpm]), 
    i.e. obtained by the formula `floor((i+2)^2 / 4)` for some ``i\\ge 0``. 
    Otherwise, the result will be incorrect.

!!! note "Remark 3"
    The active (diagonal) Padé implementation uses the closed-form formula
    ``2(U\\_e-U\\_o)^{-1}U\\_o + I`` for `:diagonalcheap`.

# References 
> [^hf19_mpexpm] M. Fasi and N. J. Higham, An Arbitrary Precision Scaling and Squaring Algorithm for the Matrix Exponential
> SIAM J. Matrix Anal. Appl., Vol. 40.4 (2019), pp.1233-1256.
> [doi: 10.1137/18M1228876](https://doi.org/10.1137/18M1228876)
"""
function eval_pade! end



########## Ricalcolo delle diagonali per lo squaring triangolare superiore ##########

"""
    recompute_diagonals!(T, Y)

Recomputes the main diagonal and first upper diagonal elements of 
``Y ≈ e^T`` using those of ``T``. 
Namely, the aforementioned elements of ``Y`` are replaced with 
quantities computed using exact formulas, computed on the original elements
(those of ``T``).

!!! note
    This function overwrites `Y`.
"""
function recompute_diagonals!(
    T::AbstractMatrix, 
    Y::AbstractMatrix
) 
    n = LinearAlgebra.checksquare(T)
    i = 1
    while i ≤ n
        # invariant: [i,i] is the top-left corner of a block 
        # (be it 1x1, 2x2, or two successive 1x1s (aka triangular 2x2))
        if (i+1 == n) || (i ≤ n-2 && T[i+2,i+1] == 0)
            # we're in a 2x2 block, eventually triangular: 
            # in such a block, T[i+2,i+1] is always 0    
            
            if T[i+1,i] == 0   # the block is triangular
                Y[i:i+1,i:i+1] = expm2by2_tri(T[i:i+1,i:i+1])
            else               # the block is full 
                Y[i:i+1,i:i+1] = expm2by2_full(T[i:i+1,i:i+1])
            end
            i += 2      # exit this block and actually go to the next one
        else 
            # we're in a 1x1 block followed by a 2x2 full block (or at the boundary)
            Y[i,i] = exp(T[i,i])
            i += 1
        end
    end
    return nothing
end


sinch(x::Real) = x == 0 ? 1 : sinh(x) / x

function sinch(z::Complex)
    if real(z) == 0 
        # Remark: we return a real result!
        return z == 0 ? 1 : imag(sinh(z)) / imag(z)
    else 
        return z == 0 ? 1 : sinh(z) / z
    end
end


"""
    Y = expm2by2_full(B)

Computes ``Y = e^B``, the exponential of a full ``2\\times 2`` 
block `B`, using formula (2.2) from [^alhi09n].

> [^alhi09n] N. J. Higham and A. H. Al-Mohy, A New Scaling and Squaring Algorithm for the Matrix Exponential
> SIAM J. Matrix Anal. Appl., Vol 31.3 (2010), pp.970-989
> [doi: 10.1137/09074721X](https://doi.org/10.1137/09074721X)
"""
function expm2by2_full(B)
    Y = similar(B)
    b11, b21, b12, b22 = B[:]
    b11mb22 = b11 - b22
    μsq = (b11mb22)^2 + 4*b12*b21
    if μsq < 0 
        μsq = Complex(μsq)
    end
    δ = sqrt(μsq)/2     # μ/2 in the formula
    exp_apd2 = exp((b11+b22)/2)
    # the result is real, but it may be of Complex type
    coshδ = real(cosh(δ))
    sinchδ = real(sinch(δ))

    Y[1,1] = exp_apd2 * (coshδ + (b11mb22)/2 * sinchδ)
    Y[2,1] = exp_apd2 * b21 * sinchδ
    Y[1,2] = exp_apd2 * b12 * sinchδ
    Y[2,2] = exp_apd2 * (coshδ + (-b11mb22)/2 * sinchδ)   
    return Y
end


"""
    expm2by2_tri(T)

Computes ``Y = e^T``, the exponential of a upper triangular ``2\\times 2``
block `T`, using formula (10.42) from [^Higham].

> [^Higham] Higham, N. J. Functions of Matrices, SIAM, 2008
> [doi:10.1137/1.9780898717778](https://doi.org/10.1137/1.9780898717778)
"""
function expm2by2_tri(M::AbstractMatrix{T}) where {T}
    Y = zero(M)
    M₁, M₂ = diag(M)

    Y[1,1] = exp(M₁)
    Y[2,2] = exp(M₂)

    M₁ += M₂        # M₁ ← M[1,1] + M[2,2]
    M₂ -= M[1,1]    # M₂ ← M[2,2] - M[1,1]

    exp_arg = M₁ / 2
    sinch_arg = M₂ / 2

    if max(real(exp_arg), abs(sinch_arg)) < log(floatmax(real(T)))    # guard against overflow
        Y[1,2] = M[1,2] * exp(exp_arg) * sinch(sinch_arg)
    else
        # Numerical cancellation if M[2,2] ≈ M[1,1] 
        # we use divided differences in this case
        Y[1,2] = M[1,2] * (Y[2,2] - Y[1,1]) / M₂
    end
    return Y
end



############ Evaluation of the upper bound on the forward error ############

function scalar_error_tayl!(
    S::AandPowsStruct,
    x,
    m::Integer,
    s::Real,
    extra_precision::Bool,
    factorials::FactorialsStruct,
    ψ,
    compute_ψ::Bool = true
)
    ν = ceil(typeof(m), √m)     # the "batch size" of P.-S.
    for i in length(S.powers)+1:ν+1
        push!(S.powers, S.powers[i-1] * S.A)    # Add powers
    end

    old_prec = precision(BigFloat)
    intrnl_prec = extra_precision ? floor(Int64, old_prec * p_factor) : old_prec

    δ = setprecision(BigFloat, intrnl_prec) do 
        xb = big(x) 
        mb = big(m)
        tₘ = sum(xb.^(0:mb) ./ factorials(0:mb))  
        δ = abs(tₘ - exp(xb))  # error bound |tₘ(α) - exp(α)|
        δ
    end
    κ_A = 1

    # Stop recomputing ψ once successive approximations agree to √eps(Float64)
    if compute_ψ
        lm1 = length(S.powers) - 1
        numerators = [2.0^(-s*k) for k in 0:lm1]
        factorials_double = factorials(0:lm1, return_type=Float64)
        coeffs = numerators ./ factorials_double
        approx = sum(coeffs .* S.powers)    
        ψ_new = opnorm(approx, 1)   
                                    
        if ψ_new ≈ ψ
            compute_ψ = false 
        end
    else 
        ψ_new = ψ
    end
                            
    return δ, ψ_new, κ_A, compute_ψ
end


function scalar_error_pade!(
    S::AandPowsStruct,
    x,
    m::Integer,
    s::Real,
    extra_precision::Bool
)
    n = LinearAlgebra.checksquare(S.A)
    T = eltype(S.A)
    
    old_prec = precision(BigFloat)
    intrnl_prec = extra_precision ? floor(Int64, old_prec * p_factor) : old_prec

    δ, ψ, κ_A = setprecision(BigFloat, intrnl_prec) do
        xb = big(x) 
        mb = big(m)
        ranges = 0:mb
        
        c_num = (factorial(mb)/factorial(2mb)) ./
                (factorial.(mb .- ranges) .* factorial.(ranges)) .*
                factorial.(2mb .- ranges)
        c_num[1] = big(1.)
        c_den = ((-1).^ranges) .* c_num    

        xx = similar(c_num)
        xx[1] = one(xb)
        for k in 2:m+1
            xx[k] = xb * xx[k-1]
        end

        out_T = T <: Real ? Float64 : ComplexF64
        Uₑ = polyvalm_ps!(S.powers, s, c_num[1:2:end], outputclass=out_T)
        if m ≥ 1 
            Uₒ = polyvalm_ps!(S.powers, s, c_num[2:2:end], outputclass=out_T)
            Uₒ = (convert(Matrix{out_T}, S.A) / 2^s) * Uₒ
        else 
            Uₒ = zero(S.A)
        end

        Qₘ = Uₑ - Uₒ

        F = lu(Qₘ) 
        opinv = InverseMap(F)

        # ‖ qₘ(2^(-s)A)^(-1) ‖₁ ≡ opnorm(inv(Qₘ), 1)
        η = opnorm(F.U \ inv(F.L) * F.P, 1)
        
        δ = η * abs(exp(xb) * dot(c_den, xx) - dot(c_num, xx))
        ψ = opnorm1est(2 * opinv * Uₒ + I(n))   # ‖ exp(2^(-s)A) ‖ = ‖ rₘ(2^(-s)A) ‖₁ ≈ ‖ exp(2^(-s)A) ‖₁
        κ_A = η * opnorm1est(LinearMap(Qₘ)) 

        δ, ψ, κ_A
    end                                    

    return δ, ψ, κ_A
end



"""
"""
function eval_error!(
    S::AandPowsStruct,
    x,
    m::Integer,
    s::Real,
    extra_precision::Bool,
    factorials::Union{FactorialsStruct,Nothing}=nothing,
    ψ_in=nothing,
    compute_ψ::Bool=true
)
    if S.use_taylor
        isnothing(factorials) && throw(ArgumentError("`factorials` is required when `S.use_taylor` is true."))
        isnothing(ψ_in) && throw(ArgumentError("`ψ_in` is required when `S.use_taylor` is true."))
        return scalar_error_tayl!(S, x, m, s, extra_precision, 
                                 factorials, ψ_in, compute_ψ)
    else 
        δ, ψ, κ_A = scalar_error_pade!(S, x, m, s, extra_precision)
        return δ, ψ, κ_A, true
    end
end



############ Computation of α_min ############

function _alpha!(
    alpha_vec::AbstractVector, 
    S::AandPowsStruct,
    k, m, s
)
    d = fld(1 + sqrt(4*(m+k) + 5), 2)  
    d = Int64(d)

    if alpha_vec[d+1] == 0 
        if alpha_vec[d] == 0 
            alpha_vec[d] = normest1(d, S)^(1/d)
        end
        # Proviamo a evitare il calcolo di ‖A^(d+1)‖
        low  = findfirst(!iszero, alpha_vec) # lowest index of a nonzero α
        high = findlast(!iszero, alpha_vec)  # highest index of a nonzero α
        bin_counter = false
        found_upper_bound = false
        while low < high
            if low + high == d + 1
                dp1od = (d+1)/d     # com'è sul MATLAB originale non mi torna
                if (alpha_vec[d])^dp1od > alpha_vec[low]*alpha_vec[high]
                    #print("Upper bound found! low = $(low), high=$(high)\n")
                    found_upper_bound = true
                    break
                end
            end
            if bin_counter
                low += 1
            else 
                high -= 1
            end
            bin_counter = !bin_counter
        end
        if found_upper_bound
            return alpha_vec[d] / 2^s
        else 
            #print("Upper bound not found...\n")
            alpha_vec[d+1] == 0 || error("Uhm qualcosa non torna qua…\n")
            alpha_vec[d+1] = normest1(d+1, S)^(1/(d+1))
        end
    end
    # oss: at one point, alpha_vec[d+1] was 0, then it has been computed.
    #      When this happened, also alpha_vec[d] had been computed.
    α_min = maximum(alpha_vec[d:d+1])   
    α_min /= 2^s
    return α_min
end


function _alpha!(
    alpha_dict::Dict,
    S::AandPowsStruct,
    k, m, s
)
    d = fld(1 + sqrt(4*(m+k) + 5), 2)   # that is d^{[k/m]}
    d = Int64(d)

    α_d = get!(alpha_dict, d) do 
        normest1(d, S)^(1/d)
    end
    
    if !haskey(alpha_dict, d+1)
        # Try to avoid computing ‖A^(d+1)‖
        found_upper_bound = false
        dp1od = (d+1) / d 
        for low_key in keys(alpha_dict)
            high_key = d + 1 - low_key
            if haskey(alpha_dict, high_key) && low_key < high_key
                if α_d^dp1od > alpha_dict[low_key]*alpha_dict[high_key]
                    found_upper_bound = true 
                    break 
                end
            end
        end
        if found_upper_bound 
            return α_d / 2^s 
        else 
            alpha_dict[d+1] = normest1(d+1, S)^(1/(d+1))            
        end        
    end
    # At one point, alpha_vec[d+1] was 0, then it has been computed.
    # When this happened, also alpha_vec[d] had been computed.
    α_min = max(alpha_dict[d+1], alpha_dict[d])
    α_min /= 2^s
    return α_min
end


function alpha! end 

function alpha!(alphaVecOrDict, S::AandPowsStruct, m, s)
    #real(eltype(S.A)) == BigFloat && @warn "Please don't use alpha! with arbitrary precision data"

    q = S.use_taylor ? 0 : m
    return _alpha!(alphaVecOrDict, S, m, q, s)
end



############ Evaluation (better: estimate) of || A^d ||₁ ############

# Function that returns the action X -> AᵈX, using only elements in `Apows`
function evalPowVecDiag end

# Taylor implementation
function evalPowVecDiag(
    d::Integer,
    S::AandPowsStruct,
    ::Val{true},
)    
    # determine the type of the matrices' elements
    T = promote_type(eltype.(S.powers)...)
    T_low = T <: Complex ? ComplexF64 : Float64

    n = LinearAlgebra.checksquare(S.A)

    Ad_action = function (X::AbstractVecOrMat)
        p = d
        l = min(length(S.powers), p+1)
        while p > 0 
            for _ in 1:fld(p, l-1)
                A_l_double = convert(Matrix{T_low}, S.powers[l])
                X = A_l_double * X
            end
            p = mod(p, l-1)         
            l = min(l-1, p+1)
        end 
        return X
    end
    Adp_action = function (X::AbstractVecOrMat)
        p = d
        l = min(length(S.powers), p+1)
        while p > 0 
            for _ in 1:fld(p, l-1)
                A_l_double = convert(Matrix{T_low}, S.powers[l])
                X = A_l_double' * X
            end
            p = mod(p, l-1)         
            l = min(l-1, p+1)
        end 
        return X
    end

    kwargs = (
        issymmetric = issymmetric(S.A),
        ishermitian = ishermitian(S.A),
        isposdef = isposdef(S.A)
    )
    return LinearMap{T_low}(Ad_action, Adp_action, n; kwargs...)
end

# Padé implementation
function evalPowVecDiag(
    d::Integer, 
    S::AandPowsStruct,
    ::Val{false}
)    
    # determine the type of the matrices' elements
    T = promote_type(eltype.(S.powers)...)
    T_low = T <: Complex ? ComplexF64 : Float64
    
    n = LinearAlgebra.checksquare(S.A)


    Ad_action = function (X::AbstractVecOrMat)
        p = d
        l = length(S.powers)
        while p > 1 && l > 1
            for _ in 1:fld(p, 2*(l-1))
                A_l_double = convert(Matrix{T_low}, S.powers[l])
                X = A_l_double * X
            end 
            p = mod(p, 2*(l-1))         # d -= ⌊d/2(l-1)⌋ * 2(l-1)
            l = min(l-1, fld(p,2)+1)
        end
        if p == 1 
            # extra multiplication by A, in case `d` was odd
            A_double = A_l_double = convert(Matrix{T_low}, S.A)
            X = A_double * X
        end    
        return X
    end
    Adp_action = function (X::AbstractVecOrMat)
        p = d
        l = length(S.powers)
        while p > 1 && l > 1
            for _ in 1:fld(p, 2*(l-1))
                A_l_double = convert(Matrix{T_low}, S.powers[l])
                X = A_l_double' * X
            end 
            p = mod(p, 2*(l-1))         # d -= ⌊d/2(l-1)⌋ * 2(l-1)
            l = min(l-1, fld(p,2)+1)
        end
        if p == 1 
            # extra multiplication by A', in case `d` was odd
            A_double = convert(Matrix{T_low}, S.A)
            X = A_double' * X
        end    
        return X
    end

    kwargs = (
        issymmetric = issymmetric(S.A),
        ishermitian = ishermitian(S.A),
        isposdef = isposdef(S.A)
    )
    return LinearMap{T_low}(Ad_action, Adp_action, n; kwargs...)
end

evalPowVecDiag(d::Integer, S::AandPowsStruct) = evalPowVecDiag(d, S, Val(S.use_taylor))


"""
    normest1(d, S)

Computes an estimate of the 1-norm of ``A^d``. If ``A`` is a real nonnegative 
matrix, the estimate is exact, otherwise the elements in `S.powers` are used.

The implementation of the action of ``A^d`` is the same as the 
EvalPowVecDiag function in [^hf19_mpexpm]. Precisely, a `LinearMap` object
that implements the action is created, then the norm is estimated using 
the `opnorm1est` function from `MatrixEquations.jl`.

!!! note
    The norm is computed in double precision. There are `convert`s 
    inside the linear operator whose norm is computed.

# Arguments
- `d::Integer`: the power of `A`
- `S::AandPowsStruct`: a struct with the fields `use_taylor`,
                       `A` (the matrix) and `powers` 
    - `A::AbstractMatrix`: the matrix ``A``
    - `use_taylor::Bool`: whether the chosen approximant is the Taylor one or not.                          
    - `powers::AbstractVector`: a vector with the powers of ``A``. 
                               It is assumed that `powers` is ``[I, A, \\dots, A^l]`` 
                               in the former case, and that it's ``[I, A^2, \\dots, A^(2l)]``
                               in the latter case.

# References 
> [^hf19_mpexpm] N. J. Higham and M. Fasi, An Arbitrary Precision Scaling and Squaring Algorithm for the Matrix Exponential
> SIAM J. Matrix Anal. Appl., Vol. 40.4 (2019), pp.1233-1256.
> [doi: 10.1137/18M1228876](https://doi.org/10.1137/18M1228876)
"""
function normest1(d::Integer, S::AandPowsStruct)
    length(S.powers) > 1 || throw(ArgumentError(lazy"Supply at least the 0th and 1st power of the matrix"))
    n = LinearAlgebra.checksquare(S.A)

    if S.A == abs.(S.A)
        e = ones(n)
        for _ in 1:d
            e = (S.A)' * e
        end
        γ_d = norm(e, Inf)
    else
        Ad_action = evalPowVecDiag(d, S)
        γ_d = MatrixEquations.opnorm1est(Ad_action)    
    end

    return γ_d
end



############ Optimal degrees ############

function opt_degs end

function opt_degs(::Val{:taylor}, max_deg::Integer=500)
    # degs[i] = ⌊(i+2)²/4⌋
    degs = [1,    2,    4,    6,    9,   12,   16,   20,   25,
            30,   36,   42,   49,   56,   64,   72,   81,   90,  100,
            110,  121,  132,  144,  156,  169,  182,  196,  210,  225,
            240,  256,  272,  289,  306,  324,  342,  361,  380,  400,
            420,  441,  462,  484,  506,  529,  552,  576,  600,  625,
            650,  676,  702,  729,  756,  784,  812,  841,  870,  900,
            930,  961,  992, 1024, 1056, 1089, 1122, 1156, 1190, 1225,
            1260, 1296, 1332, 1369, 1406, 1444, 1482, 1521, 1560, 1600,
            1640, 1681, 1722, 1764, 1806, 1849, 1892, 1936, 1980, 2025,
            2070, 2116, 2162, 2209, 2256, 2304, 2352, 2401, 2450, 2500]
    filter!(x -> x < max_deg, degs) 
end

function opt_degs(::Val{:diagonalcheap}, max_deg::Integer=500)
    # degs[i] = 2⋅⌈(i-1)/4⌉⋅(i-1-2⌊(i-2)/4⌋) + 1
    degs = [1,    2,    3,    5,    7,    9,    13,   17,   21,
            25,   31,   37,   43,   49,   57,   65,   73,   81,   91,
            101,  111,  121,  133,  145,  157,  169,  183,  197,  211,
            225,  241,  257,  273,  289,  307,  325,  343,  361,  381,
            401,  421,  441,  463,  485,  507,  529,  553,  577,  601,
            625,  651,  677,  703,  729,  757,  785,  813,  841,  871,
            901,  931,  961,  993,  1025, 1057, 1089, 1123, 1157, 1191,
            1225, 1261, 1297, 1333, 1369, 1407, 1445, 1483, 1521, 1561,
            1601, 1641, 1681, 1723, 1765, 1807, 1849, 1893, 1937, 1981,
            2025, 2071, 2117, 2163, 2209, 2257, 2305, 2353, 2401, 2451,
            2501]
    filter!(x -> x < max_deg, degs)
end



############ check that A is in Schur form ############

function isschur end

function isschur(
    A::AbstractMatrix{T},
    _atol::Real = eps(real(zero(T)))
) where {T<:Complex}
    istriu(A)
end

function isschur(
    A::AbstractMatrix{T}, 
    atol::Real = eps(real(float(T)))
) where {T<:Real}
    n = LinearAlgebra.checksquare(A)
    istriu(A) && return true 

    # return false if there are nonzero elements below the 1st subdiagonal
    tril(A, -2) == zero(A) || return false

    # check 2×2 block structure 
    i = 2
    while i ≤ n
        # invariant: [i-1,i-1] is the top-left corner of a block (1x1 or 2x2)
        if abs(A[i,i-1]) > atol
            # must be a 2×2 block (skip next index)
            if i < n && A[i+1,i] != 0   #abs(A[i+1,i] > atol)
                # there's not a zero below the bottom-right element 
                return false  
            end

            # check that the eigenvalues are complex conjugate 
            a, c, b, d = A[i-1:i,i-1:i][:]  # A[i-1:i,i-1:i] = [a b; c d]
            
            Δ = (a-d)^2 + 4*b*c
            if Δ ≥ -atol
                # the eigenvalues are complex conjugate ⟺ Δ < 0
                return false 
            end

            i += 2
        else
            i += 1
        end
    end

    return true
end

@inline update_epsilon(ϵ, factor, ::Val{false}) = ϵ * factor   
@inline update_epsilon(ϵ, factor, ::Val{true}) = ϵ

@inline my_eps(::Type{T}) where {T} = eps(real(big(float(T))))
