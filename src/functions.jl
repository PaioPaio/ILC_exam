function createsystem()
    @parameters m κ σ  # mass, stiffness, damping
    ps = [m;κ;σ]
    @variables x(t)[1:N] v(t)[1:N] u(t)

    ### System's Equations
    eqs1 = [D(x[i]) ~ v[i] for i in eachindex(x)]
    eqs2 = [D(v[1]) ~ -2κ / m * x[1] + κ / m * x[2]]
    eqs3 = [
        D(v[i]) ~ κ / m * x[i-1] - 2κ / m * x[i] + κ / m * x[i+1] + σ / m * v[i-1] - 2σ / m * v[i] + σ / m * v[i+1] for i in 2:length(x)-1
    ]
    eqs4 = [D(v[end]) ~ κ / m * x[end-1] - κ / m * x[end] + σ / m * v[end-1] - σ / m * v[end] + u / m]
    append!(eqs1, eqs2, eqs3, eqs4)

    return ps, ODESystem(eqs1, t)
end

# Defined just because exp does not work with Matrix{Num}
function expmapproxint(A,t;order=5) 
    T = promote_type(eltype(A), typeof(t))
    E = Matrix{T}(I,size(A)...)
    term = E

    for i in 1:order 
        term *= A*t/i 
        E += simplify(term)
    end

    return simplify(E)
end

expmapprox(A::AbstractMatrix{Num},t::Real) = expmapproxint(A,t)
expmapprox(A::AbstractMatrix{Real},t::Num) = expmapproxint(A,t)
expmapprox(A::AbstractMatrix{Num},t::Num) = expmapproxint(A,t)
expmapprox(A,t) = exp(A*t)

## Function Basis
function CAtB(A,B,C,t)
    return C * expmapprox(A,t) * B
end
CAtB(A::AbstractMatrix{T},B::AbstractMatrix{T},C::AbstractMatrix{T},t::T) where T<:AbstractFloat = C*exp(A*t)*B

function pif(A,B,C,T,i,t)
    minT = i==1 ? 0 : T[i-1]
    maxT = T[i]
    return ifelse(!((t>maxT)&(t<minT)), CAtB(A,B,C,maxT-t)', zeros(size(B,2),size(C,1)))
    #return ifelse(!((t>T[i])&(t<T[i-1])), B' * (exponential(A'*(T[i]-t))) * C', zeros(size(B,2),size(C,1)))
end
