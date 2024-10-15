using OrdinaryDiffEq, LinearAlgebra, Integrals, GLMakie

includet("../src/fun2.jl")

# System parameters
N = 5
T = [8.0, 18.0, 20.0]
o = length(T)
m, κ, σ = 1.0, 2.0, 1.0
pp = (m,κ,σ)
Niter = 30

# Define system matrices

A, B = system_matrices(N, m, κ, σ)
C = float(I(2N))  # y = [x; v], so C is identity

##
# Desired output
x₁ = cumsum(ones(N))./N
y₁ = vcat(x₁, zeros(N))
ywant = hcat(y₁, zeros(2N), zeros(2N))

# Initialize storage for results
xⁱ = zeros(N, o, Niter)                         
vⁱ = zeros(N, o, Niter)
uⁱ = zeros(1, o, Niter)
αrray = zeros(2N, o, Niter)
solarray = Vector{Any}(undef, Niter)

##  
# Compute H matrix

H = compute_H(A, B, C, T)
L = inv(H'H + I(2N*o)*1e-2) * H'

## ILC loop
u0 = vcat(zeros(N), zeros(N), 0.0)
tspan = (0.0, T[end])
for i in 1:Niter
    # Set up and solve ODE
    p = (A, B, C, αrray[:,:,i], T,utotn)
    prob = ODEProblem(mass_spring_damper2!, u0, tspan, p)
    sol = solve(prob,Tsit5(),abstol= 1e-9)
    solarray[i] = sol

    # Extract and store results
    for (j, t) in enumerate(T)
        xⁱ[:, j, i] = sol(t)[1:N]
        vⁱ[:, j, i] = sol(t)[N+1:2N]
        #uⁱ[:, j, i] = sol(t)[end:end]
    end

    # Update α for next iteration
    if i < Niter
        y_current = vcat(xⁱ[:,:,i], vⁱ[:,:,i])
        αrray[:,:,i+1] = αrray[:,:,i] + reshape(L * (ywant[:] - y_current[:]), 2N, o)
    end
end

##
# Time vector from the solution (full time points)
t_vals = 0.0:0.05:T[end] 
x_vals = [solarray[end](t)[1:N] for t in t_vals]  # Position data
v_vals = [solarray[end](t)[N+1:2N] for t in t_vals]  # Velocity data 
u_vals = [utotn(A, B, C, αrray[:,:,end], t, T) for t in t_vals]       
#u_vals = [solarray[end](t)[end] for t in t_vals]  # Control input over time

# Creating the figure
f = Figure(resolution = (1200, 1000))

# 1. Position plot
ax1 = Axis(f[1, 1], title="Positions", xlabel="Time", ylabel="Positions x₁ to xₙ")
for i in 1:N
    lines!(ax1, t_vals, getindex.(x_vals, i), label="x$i")  # Positions
    scatter!(ax1, T, ywant[i,:], color=:black, marker=:x, markersize=8)  # Desired positions as x markers
end
axislegend(ax1)

# 2. Velocity plot
ax2 = Axis(f[2, 1], title="Velocities", xlabel="Time", ylabel="Velocities v₁ to vₙ")
for i in 1:N
    lines!(ax2, t_vals, getindex.(v_vals, i), label="v$i")  # Velocities
    scatter!(ax2, T, ywant[N+i,:], color=:black, marker=:x, markersize=8)  # Desired velocities as x markers
end
axislegend(ax2)
#
base_vals = zeros(o,2N,length(t_vals))
# 3. Base functions plot using upi
ax3 = Axis(f[3, 1], title="Base Functions π", xlabel="Time", ylabel="Base Functions")
for i in 1:o
    base_vals[i,:,:] .= hcat([upi(A, B, C, t, T, i) for t in t_vals]...)
    for j in 1:2N# Using upi to compute the base functions
        lines!(ax3, t_vals, base_vals[i,j,:])
    end
end
#
# 4. Control Action plot (input uₙ)
ax4 = Axis(f[4, 1], title="Control Action uₙ", xlabel="Time", ylabel="Control Action")
lines!(ax4, t_vals, u_vals, color=:black)

# Display the figure
display(f)