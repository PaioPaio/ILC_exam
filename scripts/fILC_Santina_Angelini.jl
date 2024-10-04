using DrWatson, ModelingToolkit, GLMakie
using ModelingToolkit: t_nounits as t, D_nounits as D

N = 30
@parameters m κ σ  # mass, stiffness, damping
@variables x(t)[1:N] v(t)[1:N] u(t)
##
eqs1 = [D(x[i]) ~ v[i] for i in eachindex(x)]
eqs2 = [D(v[1]) ~ -2κ / m * x[1] + κ / m * x[2]]
eqs3 = [
    D(v[i]) ~ κ / m * x[i-1] - 2κ / m * x[i] + κ / m * x[i+1] + σ / m * v[i-1] - 2σ / m * v[i] + σ / m * v[i+1] for i in 2:N-1
]
eqs4 = [D(v[end]) ~ κ / m * x[end-1] - κ / m * x[end] + σ / m * v[end-1] - σ / m * v[end] + u / m]
append!(eqs1, eqs2, eqs3, eqs4)
##
@named sys = ODESystem(eqs1, t)
f, x_sym, ps = ModelingToolkit.generate_control_function(sys, [u], simplify=true)

###
p = [1, 2, 1] # mass, stiffness, damping
x₀ = zeros(N)
v₀ = zeros(N)
x₁ = cumsum(ones(N) ./ N)
x₂ = zeros(N)
x₃ = zeros(N)

###
