"""
    AandPowsStruct(A, Apows, use_taylor)

A struct to hold the data used by the `exp_mp` internal functions. 
In particular, it contains the powers of the matrix

# Fields 
- `A::AbstractMatrix`: The base matrix
- `powers::AbstractVector{<:AbstractMatrix}`: The powers of `A` or `A^2` (depending on `use_taylor`)
- `use_taylor::Bool`: Whether the struct is used by Taylor-related functions or not 

# Notes 
Besides the type specifications of the struct's record, there is no check. 
In particular, it's not checked that `A` is square, that `powers` and `use_taylor`
are consistent, that the matrices in `powers` are the right ones. 
**This saves a lot of overhead, but the responsibility for a correct use is on the 
implementation**
"""
struct AandPowsStruct
    A::AbstractMatrix
    powers::AbstractVector{<:AbstractMatrix}
    use_taylor::Bool
end



"""
    FactorialsStruct(f_vec)
    FactorialsStruct(n::Integer)
    FactorialsStruct()

A callable struct used to cache the values of factorials for repeated access.

# Fields
- `f_vec::Vector{BigInt}`: The cached factorials, with `f_vec[k+1] = k!` for `k >= 0`.

"""
struct FactorialsStruct
    f_vec::Vector{BigInt}
    FactorialsStruct(n::Integer) = 
        begin
            n < 0 && throw(DomainError(lazy"FactorialsStruct expects n to be positive"))
            v = [factorial(big(k)) for k=0:n]
            new(v)
        end
    FactorialsStruct() = new([big(1)])
end


function (f::FactorialsStruct)(
    n::Integer; 
    return_type::DataType=BigInt
)
    n < 0 && throw(DomainError(lazy"FactorialsStruct call expects n to be positive"))
    l = length(f.f_vec)
    if n+1 > l
        nuovi_f = l:n .|> big .|> factorial
        append!(f.f_vec, nuovi_f)
    end
    return return_type(f.f_vec[n+1])
end


function (f::FactorialsStruct)(
    r::AbstractRange{<:Integer}; 
    return_type::DataType=BigInt
)
    fst, lst = first(r), last(r)
    stp = step(r)

    (fst < 0 || lst < 0) && throw(DomainError(lazy"FactorialsStruct call expect ranges that span positive integers"))
    isempty(r) && return return_type[]
    
    l = length(f.f_vec)
    mx = max(fst, lst)
    if mx+1 > l
        nuovi_f = l:mx .|> big .|> factorial
        append!(f.f_vec, nuovi_f)
    end
    return return_type.(f.f_vec[fst+1:stp:lst+1])
end
