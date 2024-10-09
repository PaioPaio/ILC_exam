using DrWatson, ModelingToolkit, OrdinaryDiffEq, Symbolics, LinearAlgebra, Integrals
using ModelingToolkit: t_nounits as t, D_nounits as D   #you can define these yourself, but included

includet("../src/functions.jl")

## System Parameters
N = 5
o = 3
x₀ = zeros(N)
v₀ = zeros(N)
ps = @parameters m κ σ α[1:2N,1:o] # mass, stiffness, damping
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
p = [1, 2, 1, [zeros(2*N,o)]] # mass, stiffness, damping
Niter = 100

## System's Equations
eqs = [ [D(x[i]) ~ v[i] for i in eachindex(x)]          #dotx=v
        [D(v[1]) ~ -2κ / m * x[1] + κ / m * x[2]]       #dotdotx = stuff
        [D(v[i]) ~ κ / m * x[i-1] - 2κ / m * x[i] + κ / m * x[i+1] + σ / m * v[i-1] - 2σ / m * v[i] + σ / m * v[i+1] for i in 2:length(x)-1]
        [D(v[end]) ~ κ / m * x[end-1] - κ / m * x[end] + σ / m * v[end-1] - σ / m * v[end] + u / m]
        [y[i] ~ x[i] for i in 1:N]                      #y=x
        [y[i] ~ v[i-N] for i in N+1:2N]]                
## Define System without input policy and obtain jacobian/state matrices (just because I'm lazy)
@named sysnoinput = ODESystem(eqs, t)
J = calculate_jacobian(sysnoinput)

A = J[1:2N, 1:2N]
B = J[1:2N, 2N+1]
C = J[2N+1:4N,1:2N]

An = Float64.(substitute(A, Dict(ps.=> p)))
Bn = Float64.(substitute(B, Dict(ps.=> p)))
Cn = Float64.(substitute(C, Dict(ps.=> p)))     

##  
# Parametrized Inputs (basis defined in src file)

function utotn(α::AbstractVector{<:Num}, s::Num) ::Num
    return uiter(An, Bn, Cn, T, s, α)
end

@register_symbolic utotn(α,s)
##              
eqsinput = [u ~ utotn(α,t)]
## fILC
# essentially @named but better as it simplifies stuff along the way if it can. Before we could not use it as u wasnt defined
#@named sys = ODESystem(push!(eqs,eqsinput),t)

eqsnew = vcat(eqs,eqsinput)
@mtkbuild sys = ODESystem(eqsnew,t)

##
xⁱ = zeros(N,o,Niter)
vⁱ = zeros(N,o,Niter) 
solarray = Vector{Any}(undef,Niter)
αrray = zeros(2*N,o,Niter)

method = N<=5 ? HCubatureJL() : VEGASMC()

## Build H
Twith0 = [0.0 T[:]...]
ftoint1 = [i>j ? ((t,p)-> Cn*exp(An*(Twith0[i]-t))*Bn*Bn'* exp(An'*(Twith0[j]-t))*Cn') : zeros(2N,2N)
            for i in 2:o+1, j in 1:o]

# subdivide in multiple integrals so that it's easier (pis have bounded domain)
H = zeros(2N*o,2N*o)
for i in 1:o, j in 1:o 
    blockrangei = (i-1)*2N+1:i*2N
    blockrangej = (j-1)*2N+1:j*2N

    H[blockrangei,blockrangej] .= i>=j ? solve(IntegralProblem(ftoint1[i,j],(Twith0[j],Twith0[j+1])),method).u[:,:] : ftoint1[i,j]
end

L = inv(H'H+I(2N*o)*1e-2)*H'

##

for i in 1:Niter
    # Instantiate the symbolic System
    systosolve = ODEProblem(sys, #system to solve
                            [x=>x₀,v=>v₀], #control with α to iterate over
                            (T[0], T[end]), #time span
                            ps.=>p) #parameters (includes different α)
    sol = solve(systosolve, Tsit5())
    
    H = vcat([solve(Hprob[i],method) for i in eachindex(Hprob)]...)
    solarray[i] = sol
    xⁱ[:,:,i] = sol[x](T)
    vⁱ[:,:,i] = sol[v](T)
    αrray[:,:,i+1] = αrray[:,:,i] + reshape(L*([ywant - hcat(xⁱ[:,:,i],yⁱ[:,:,i])[:]]),2N,o)    
    p[end] = αrray[:,:,i+1]
end