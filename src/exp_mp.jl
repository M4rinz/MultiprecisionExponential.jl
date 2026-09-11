const ASSRT_STRNG = """
    Taylor is being used but epsilon^(-1/8) ≤ 1. 
    Maybe you're being too lenient with the tolerance on the error bound (i.e. epsilon is too big)?
    """



"""
    exp_mp(A, kwargs...)

Computes the exponential of ``A`` in arbitrary precision, using the algorithm described in [^hf19_mpexpm] 

# Arguments 
- `A::AbstractMatrix`: The matrix whose exponential we want to compute 

## Keyword arguments 
- `working_precision::Integer`: The working precision of the algorithm.
    Defaults to the precision of the type of the elements of ``A``
- `epsilon`: tolerance on the error bound used in the parameter selection phase. 
    If set to `nothing`, `epsilon` will equal the unit roundoff of the type of the elements of ``A``.
    Defaults to `nothing`
- `maxscaling::Integer`: the maximum allowed value for the sclaing parameter ``s``. Defaults to ``100``
- `maxdegree::Integer`: the maximum allowed approximant degree ``m``. Defaults to ``500``
- `algorithm::Symbol`: the algorithmic variant to use.
    It must be one of `:transfree`, `:realschur` or `:complexschur`. 
    The `:realschur` and `:complexschur` variants compute a (real and complex, respectively) Schur form ``A=QXQ^H``;
    the algorithm is then executed on ``X``, thereby computing an approximation ``R`` to ``\\exp(X)``. 
    Finally, ``QRQ^H`` is returned. The `:transfree` variant computes no Schur form.
    Note that the distinction between `:realschur` and `:complexschur` only makes sense if ``A`` is real. See remark below.
    Defaults to `:transfree`.
- `approximant::Symbol`: The approximant to use. 
    It must be either `:taylor` or `:diagonalcheap`, corresponding to Taylor and diagonal Padé 
    approximants, respectively. Defaults to `:diagonalcheap`
- `use_abs_err::Bool`: whether to use the absolute error in the stopping criterion.
    Defaults to `false`.
- `use_extra_precision::Bool`: whether to increase the precision in some internal computations 
    (such as evaluating matrix polynomials) performed when the upper bound is evaluated. 
    The local increase is by a factor of $(p_factor).
    Defaults to `true`.

!!! note "Remarks"
    - The output matrix will be of type `BigFloat` (or `Complex{BigFloat}`, if ``A`` is complex) 
    - Note that `:realschur` and `:complexschur` have the same behaviour if ``A`` is complex, 
        as `schur(A)` computes the complex Schur form by default in this case. 
        Thus, the distinction only makes sense for real matrices.

# Example
```jldoctest
julia> using LinearAlgebra, Random

julia> Random.seed!(42)

julia> A = [big(π)im 0; 0 big(π)im]
2×2 Matrix{Complex{BigFloat}}:
 0.0+3.14159im  0.0+0.0im
 0.0+0.0im      0.0+3.14159im

julia> Q, _ = qr(randn(BigFloat, 2,2));

julia> Y = exp_mp(Q*A*Q');

julia> norm(Q'*Y*Q + I(2)) ≤ 10eps(BigFloat)
true
```

# References 
> [^hf19_mpexpm] M. Fasi and N. J. Higham, An Arbitrary Precision Scaling and Squaring Algorithm for the Matrix Exponential
> SIAM J. Matrix Anal. Appl., Vol. 40.4 (2019), pp.1233-1256.
> [doi: 10.1137/18M1228876](https://doi.org/10.1137/18M1228876)
"""
function exp_mp(
    A::AbstractMatrix{T};
    working_precision::Integer = precision(real(float(T))),
    epsilon::Union{AbstractFloat,Nothing} = nothing,    # tolerance on the error bound
    maxscaling::Integer = 100,
    maxdegree::Integer = 500,
    algorithm::Symbol = :transfree,         # schur?
    approximant::Symbol = :diagonalcheap,   # chosen approximant
    use_abs_err::Bool = false,              # as stopping crit. in the main while loop
    use_extra_precision::Bool = true        # increase precision in some internal computations?
) where {T<:Number}
    n = LinearAlgebra.checksquare(A)

    approximant in VALID_APPRX || 
        throw(ArgumentError("Invalid approximant: $approximant. Valid options are $(VALID_APPRX)"))
    algorithm in VALID_ALGS || 
        throw(ArgumentError("Invalid algorithm: $algorithm. Valid options are $(VALID_ALGS)"))

    # Remark: if A is double, BigFloat's precision will match that of doubles
    old_prec = precision(BigFloat)
    try
    setprecision(BigFloat, working_precision)
    # Cast the input into (Complex) BigFloat
    A = T <: Real ? BigFloat.(A) : Complex{BigFloat}.(A)  

    if isnothing(epsilon)
        epsilon = my_eps(T)
    end
    ζ = epsilon^(-1/8)  
    use_abs_err_flag = Val(use_abs_err)

    compute_schur = false
    cmplx_schur   = true
    recompute_diag_blocks = true

    if algorithm === :transfree
        if ishermitian(A)
            d, V = eigen(A)     
            Y = V * diagm(exp.(d)) / V
            return Y
        end
        X = copy(A)
        if isschur(A) 
            recompute_diag_blocks = true
        else 
            recompute_diag_blocks = false 
        end
    elseif algorithm === :realschur
        if isschur(A)
            X = copy(A)
        else 
            compute_schur = true 
            cmplx_schur = false 
        end
    else 
        if istriu(A)    # A is already in (complex) Schur form
            X = copy(A)
        else 
            compute_schur = true
            cmplx_schur = true
        end
    end
    if compute_schur
        F = schur(A)
        if cmplx_schur
            F = Schur{Complex}(F)
        end
        X = F.T
    end

    # shift the input matrix
    μ = tr(X) / n
    useshift = true 
    positive_shift = false
    if abs(μ) > 10
        X -= μ*I(n)
        positive_shift = real(μ) ≥ 0
    else
        useshift = false
    end
   
    alpha_vec = Dict(1 => opnorm(X, 1))

    if isdiag(X)
        Y = diagm(exp.(diag(X)))
    else 
        degrees = opt_degs(Val(approximant), maxdegree)
        maxdegree = degrees[end]    # max degree `m`
        currcost = 3            # current cost (in matmuls) to 
                                # evaluate the approximant
        m = degrees[currcost]
        s = 0 

        if approximant === :taylor
            use_taylor = true
            Xpows = [I(n), X]
            factorials = FactorialsStruct(m)
        else 
            use_taylor = false
            Xpows = [I(n), X^2]
            factorials = nothing
        end
        XandP = AandPowsStruct(X, Xpows, use_taylor)

        found_degree = false 

        ψ = typemax(float(real(T)))     
        compute_ψ = true         

        # In the Taylor case, attempt using the max degree approximant, scaling as little as possible
        if alpha_vec[1] > 1e7   # if opnorm(X, 1) > 1e7
            extra_precision = false 
            if use_taylor
                α = alpha!(alpha_vec, XandP, maxdegree, s)
                δ, ψ, _, compute_ψ = eval_error!(XandP, α, maxdegree, s, extra_precision, 
                                                factorials, ψ, compute_ψ)
                while (δ > epsilon * ψ || !isfinite(ψ)) && s ≤ maxscaling
                    s += 1
                    #X /= 2
                    compute_ψ = true
                    α = alpha!(alpha_vec, XandP, maxdegree, s)
                    δ, ψ, _, compute_ψ = eval_error!(XandP, α, maxdegree, s, extra_precision,
                                                    factorials, ψ, compute_ψ)
                end
            end
        end
        extra_precision = use_extra_precision

        α = alpha!(alpha_vec, XandP, m, s)
        δ, ψ, κ_A, compute_ψ = eval_error!(XandP, α, m, s, extra_precision,
                                          factorials, ψ, compute_ψ)
        curr_ϵ = update_epsilon(epsilon, ψ, use_abs_err_flag)

        while !isfinite(ψ)
            @info "ψ = $ψ"
            currcost += 1
            s += 1
            compute_ψ = true 
            m = degrees[currcost]
            α = alpha!(alpha_vec, XandP, m, s)
            δ, ψ, κ_A, compute_ψ = eval_error!(XandP, α, m, s, extra_precision,
                                              factorials, ψ, compute_ψ)
        end


        δ_old = typemax(float(real(T)))
        
        ## main loop
        while !found_degree && m < maxdegree && s < maxscaling
            if δ ≤ curr_ϵ && κ_A < ζ
                found_degree = true 
            elseif κ_A ≥ ζ || (δ_old > 1 && abs(δ_old) < abs(δ)^2)
                assertion = κ_A ≥ ζ && use_taylor   # if Taylor is chosen, κ_A ≥ ζ should always be false
                !assertion || throw(AssertionError(ASSRT_STRNG))
                #X /= 2
                s += 1
                compute_ψ = true
            else 
                currcost += 1 
                m = degrees[currcost]
            end
            δ_old = δ
            α = alpha!(alpha_vec, XandP, m, s)
            δ, ψ, κ_A, compute_ψ = eval_error!(XandP, α, m, s, extra_precision,
                                              factorials, ψ, compute_ψ)
            curr_ϵ = update_epsilon(epsilon, ψ, use_abs_err_flag)
        end

        # The degree is now fixed. If no appropriate degree has been found, we resort to 
        # scaling as much as possible, until the upper bound on the error is < epsilon
        if !found_degree
            α = alpha!(alpha_vec, XandP, maxdegree, s)
            δ, ψ, κ_A, compute_ψ = eval_error!(XandP, α, maxdegree, s, extra_precision,
                                              factorials, ψ, compute_ψ)
            curr_ϵ = update_epsilon(epsilon, ψ, use_abs_err_flag)
            while δ ≥ curr_ϵ && s < maxscaling 
                #X /= 2
                s += 1 
                compute_ψ = true 
                α = alpha!(alpha_vec, XandP, maxdegree, s)
                δ, ψ, κ_A, compute_ψ = eval_error!(XandP, α, m, s, extra_precision,
                                                  factorials, ψ, compute_ψ)
                curr_ϵ = update_epsilon(epsilon, ψ, use_abs_err_flag)
            end
        end
        
        X /= 2^s    ## scaling 

        Y = eval_pade!(XandP, m, s, Val(approximant)) # Y = rₘ(2^(-s)X) 

        if recompute_diag_blocks
            recompute_diagonals!(X, Y)  # overwrites Y
        end
        if useshift && !positive_shift
            scal = 2.0^(-s)*μ
            Y .*= exp(scal)      # Y = e^(μ/(2^s))rₘ(2^(-s)(A-μI)) (see [exptayotf18])
            X += scal * I(n)     # X = 2^(-s)A 
        end

        ## Squaring 
        for _ in 1:s
            Y *= Y          # Y ← Y²
            if recompute_diag_blocks
                X *= 2      # X = 2^(-s+t)A if !positive_shift, else 2^(-s+t)(A-μI) 
                recompute_diagonals!(X, Y)
            end 
        end

    end # if isdiag(X)

    if compute_schur 
        Y = F.Z * Y * F.Z'
    end
    if useshift && positive_shift 
        Y .*= exp(μ)
    end
    if isreal(A)
        Y = real(Y)
    end

    return Y 
    finally #
        setprecision(BigFloat, old_prec)
    end
end



