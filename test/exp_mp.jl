using LinearAlgebra

const VALID_APPRX = (:taylor, :diagonalcheap)
const VALID_ALGS  = (:transfree, :realschur, :complexschur)


@testset "passing invalid algorithm and approximant parameter values" begin
    A = rand(2,2)

    @test_throws TypeError exp_mp(A, approximant="taylor")
    @test_throws TypeError exp_mp(A, algorithm="transfree")

    @test_throws ArgumentError exp_mp(A, approximant=:pade)
    @test_throws ArgumentError exp_mp(A, algorithm=:original)

end

@testset "random Float64 matrix" begin
    n = 5
    A = rand(n,n)
    Y_ref = exp(A)

    for appx in VALID_APPRX, alg in VALID_ALGS
        Y = exp_mp(A, algorithm=alg, approximant=appx)
        @test isapprox(Y, Y_ref)
    end
end

@testset "random ComplexF64 matrix" begin 
    n = 5
    A = rand(ComplexF64, n,n)
    Y_ref = exp(A)

    for appx in VALID_APPRX, alg in VALID_ALGS
        Y = exp_mp(A, algorithm=alg, approximant=appx)
        @test isapprox(Y, Y_ref)
    end
end

@testset "random symmetric Float64 matrix" begin
    n = 5
    A = rand(n,n)
    A = (A + A')/2
    Y_ref = exp(A)

    for appx in VALID_APPRX, alg in VALID_ALGS
        Y = exp_mp(A, algorithm=alg, approximant=appx)
        @test isapprox(Y, Y_ref)
    end
end

@testset "random symmetric ComplexF64 matrix" begin 
    n = 5
    A = rand(ComplexF64, n,n)
    A = (A + A')/2
    Y_ref = exp(A)

    for appx in VALID_APPRX, alg in VALID_ALGS
        Y = exp_mp(A, algorithm=alg, approximant=appx)
        @test isapprox(Y, Y_ref)
    end
end

@testset "Particular Schur forms" begin 
    # real matrix in real Schur form (this is normal nonSymmetric)
    A = [
         [0 1; -1 0] zeros(2,3); 
         zeros(2,2) [1 2; -2 1] zeros(2,1); 
         zeros(1,4) -1
        ]
    Y_ref = exp(A)

    for alg in VALID_ALGS
        Y = exp_mp(A, algorithm=alg)
        @test isapprox(Y, Y_ref)
    end

    # real upper triangular matrix (aka in complex Schur form)
    A = triu(randn(5,5))
    Y_ref = exp(A)

    for alg in VALID_ALGS
        Y = exp_mp(A, algorithm=alg)
        @test isapprox(Y, Y_ref)
    end

    # complex matrix in (complex) Schur form (aka upper triangular)
    R = triu(randn(ComplexF64, 5,5))
    Y_ref = exp(A)

    for alg in VALID_ALGS
        Y = exp_mp(A, algorithm=alg)
        @test isapprox(Y, Y_ref)
    end
end

@testset "Shifting" begin 
    n = 4
    A1 = 10 .* (2rand(n,n) .- 1)
    
    μ = tr(A1) / n
    A1 += (15 - μ)*I(n)
    Y_ref = exp(A1)

    for alg in VALID_ALGS, appx in VALID_APPRX
        Y = exp_mp(A1, algorithm=alg, approximant=appx)
        @test isapprox(Y, Y_ref)
    end

    A2 = copy(A1)

    μ = tr(A2) / n
    A2 += (-15 - μ)*I(n)
    Y_ref = exp(A2)

    for alg in VALID_ALGS, appx in VALID_APPRX
        Y = exp_mp(A2, algorithm=alg, approximant=appx)
        @test isapprox(Y, Y_ref)
    end    
end

@testset "Outer product at various precisions" begin 
    A = [1 2 3; 1 2 3; 1 2 3]

    function eval_Y_ref()
        divdiff = (expm1(big(3)+2+1)) / (3+2+1)
        Y_ref = I(3) + divdiff*big.(A)
        return Y_ref
    end

    for prec in [53, 256, 851, 1024],
        alg in VALID_ALGS, 
        appx in VALID_APPRX
        Y = exp_mp(A, 
            working_precision=prec, 
            algorithm=alg, 
            approximant=appx)
        setprecision(BigFloat, prec) do 
            Y_ref = eval_Y_ref()
            @test isapprox(Y, Y_ref)
        end
    end
end


