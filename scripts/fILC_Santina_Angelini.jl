using DrWatson, ModelingToolkit, OrdinaryDiffEq, GLMakie, Latexify, Symbolics, LinearAlgebra, Integrals
using ModelingToolkit: t_nounits as t, D_nounits as D   #you can define these yourself, but included

include("../src/functions.jl")

## System Parameters
N = 5
o = 3
x₀ = zeros(N)
v₀ = zeros(N)
ps = @parameters m κ σ α[1:2N*o] # mass, stiffness, damping
vars = @variables x(t)[1:N], v(t)[1:N]
inputs = @variables u(t) [input=true]
outputs = @variables y(t)[1:2N] [output=true]
## Initial and Desired Conditions
# Dimensions n=2N, l=1, m=2N, o=3
x₁ = cumsum(ones(N) ./ N)

y₁ = vcat(x₁,v₀)
y₂ = zeros(2N)
y₃ = zeros(2N)

ywant = vcat(y₁,y₂,y₃)

## Parameters
T = [8. 18. 20.] # first time is always zero
p = [1, 2, 1, [zeros(2*N*o)]] # mass, stiffness, damping
Niter = 100

## System's Equations
eqs = [ [D(x[i]) ~ v[i] for i in eachindex(x)]          #dotx=v
        [D(v[1]) ~ -2κ / m * x[1] + κ / m * x[2]]       #dotdotx = stuff
        [D(v[i]) ~ κ / m * x[i-1] - 2κ / m * x[i] + κ / m * x[i+1] + σ / m * v[i-1] - 2σ / m * v[i] + σ / m * v[i+1] for i in 2:length(x)-1]
        [D(v[end]) ~ κ / m * x[end-1] - κ / m * x[end] + σ / m * v[end-1] - σ / m * v[end] + u / m]
        [y[i] ~ x[i] for i in 1:N]                      #y=x
        [y[i] ~ v[i-N] for i in N+1:2N]]                
#append!(eqs, eqs2, eqs3, eqs4, eqsoutstate, eqsoutvel)

## Define System without input policy and obtain jacobian/state matrices (just because I'm lazy)
@named sysnoinput = ODESystem(eqs, t)
J = calculate_jacobian(sysnoinput)

A = J[1:2N, 1:2N]
B = J[1:2N, 2N+1]
C = J[2N+1:4N,1:2N]

An = Float64.(substitute(A, Dict(ps.=> p)))
Bn = Float64.(substitute(B, Dict(ps.=> p)))
Cn = Float64.(substitute(C, Dict(ps.=> p)))     

## functions

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


# Parametrized inputs
# Each elements (o=3) is a function which outputs a vector of size l×m = 1x2N
piarr = [s->pif(A,B,C,T,i,s) for i in eachindex(T)]
piarrn = [s->pif(An,Bn,Cn,T,i,s) for i in eachindex(T)]
piall(s) = [piarr[i](s) for i in eachindex(piarr)]
pitot(s)  = vcat(piall(s)...)
pitotn(s) = hcat(hcat,[piarrn[i](s) for i in eachindex(piarrn)]...)

#utot(α,s) = mapreduce(x->dot(x[1](s),x[2]),+,zip(piarr,collect(eachcol(reshape(α,2N,o)))))
#utotn(α,s) = mapreduce(x->dot(x[1](s),x[2]),+,zip(piarrn,collect(eachcol(reshape(α,2N,o)))))
#

utot(α,s) = pitot(s)*α[:]
utotn(α,s) = pitotn(s)*α
@register_symbolic utot(α,s)
##
eqsinput = u ~ utot(α,t)
##

fint(t,p) = CAtB(p[1],p[2],p[3],t)*utot(p[4],t)

## fILC
# essentially @named but better as it simplifies stuff along the way if it can. Before we could not use it as u wasnt defined
@mtkbuild sys = ODESystem(append!(eqs,eqsinput),t)

xⁱ = zeros(N,o,Niter)
vⁱ = zeros(N,o,Niter) 
solarray = Vector{Any}(undef,Niter)
αrray = zeros(2*N*o,Niter)

method = N<=5 ? HCubatureJL() : VEGASMC()

for i in 1:Niter
    # Instantiate the symbolic System
    systosolve = ODEProblem(sys, #system to solve
                            [x=>x₀,v=>v₀], #control with α to iterate over
                            (T[0], T[end]), #time span
                            ps.=>p) #parameters (includes different α)
    sol = solve(systosolve, Tsit5())
    
    Hprob = [IntegralProblem(fint,(0,T[i]),(An,Bn,Cn,p[end])) for i in 2:o+1]
    H = vcat([solve(Hprob[i],method) for i in eachindex(Hprob)]...)
    solarray[i] = sol
    xⁱ[:,:,i] = sol[x](T)
    vⁱ[:,:,i] = sol[v](T)
    αrray[:,i+1] = αrray[:,i] + pinv(H)*([ywant - hcat(xⁱ[:,:,i],yⁱ[:,:,i])[:]])    
    p[end] = αrray[:,i+1]
end