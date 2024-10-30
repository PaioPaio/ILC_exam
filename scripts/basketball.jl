using OrdinaryDiffEq, LinearAlgebra, Integrals, CairoMakie, ColorSchemes

includet("../src/functions.jl")

# System parameters
const N = 3         # in this case N can only be 3 and is the number of states (not states + velocities)
##
T = [1.0, 2.0, 3.0, 4.0, 5.5]
o = length(T)
m, w1, w2 = 0.5, 0.5, -0.9
dT = 0.5
dist = [w1; w2; -m * 9.81]
pp = (m, dist)
Niter = 30

# Define system matrices

A, B, C = basket_matrices(m)

No, Ns = size(C)        # number of observables (only the positions in this case) and of states (position + velocity)
l = size(B, 2)
##
# Desired output
ywant = zeros(3, o)
ywant[:, end] = [6.75; 0.0; 1.55]

# Initialize storage for results
xⁱ, vⁱ, uⁱ, yⁱ, αrray, solarray = initialize_sol(Ns, No, l, o, Niter)

##  
# Compute H matrix

H = compute_H_constr(A, B, C, T, dT)               # Lower block triangular as axpected
L = inv(H'H + I(N * o) * 1e-2) * H'    # Levenberg-Marquaad (or something along that way) with ϵ = 1e-2

# ILC loop
u0 = zeros(Ns)
tspan = (0.0, T[end])
for i in 1:Niter
    # Set up and solve ODE
    p = (A, B, C, αrray[:, :, i], T, dT, utotconstr, dist)
    prob = ODEProblem(LinearDynamicsDisturbed!, u0, tspan, p)
    sol = solve(prob, Tsit5(), abstol=1e-9)
    solarray[i] = sol

    # Extract and store results
    for (j, t) in enumerate(T)
        xⁱ[:, j, i] = sol(t)[1:No]
        vⁱ[:, j, i] = sol(t)[No+1:2N]
    end

    # Update α for next iteration
    if i < Niter
        yⁱ[:, :, i] = xⁱ[:, :, i]
        αrray[:, :, i+1] = αrray[:, :, i] + reshape(L * (ywant[:] - yⁱ[:, :, i][:]), No, o)
    end
end

##
# Time vector from the solution (full time points)
t_vals = 0.0:0.05:T[end]
x_vals = [solarray[end](t)[1:N] for t in t_vals]  # Position data
v_vals = [solarray[end](t)[N+1:2N] for t in t_vals]  # Velocity data 
u_vals = [utotn(A, B, C, αrray[:, :, end], t, T) for t in t_vals]

# Creating the figure
f = Figure(resolution=(1200, 1000))

plotcolors = ColorSchemes.matter

# 1. Position plot
ax1 = Axis(f[1, 1],
    title="Positions",
    xlabel="Time",
    ylabel=L"q_{i}",
    ylabelrotation=0,  # Rotates y-label 90 degrees
    ylabelsize=18,
    xticks=t_vals[1]:2.0:t_vals[end]
)
for i in 1:N
    lines!(ax1, t_vals, getindex.(x_vals, i), label=L"q_{%$i}")  # Positions
    scatter!(ax1, T, ywant[i, :], color=:black, marker=:x, markersize=8, label=L"$x$")  # Desired positions 
end
axislegend(ax1, unique=true, orientation=:horizontal)

# 2. Velocity plot
ax2 = Axis(f[2, 1],
    title="Velocities",
    xlabel="Time",
    ylabel=L"\dot{q}_{i}",
    ylabelrotation=0,  # Rotates y-label 90 degrees
    ylabelsize=18,
    xticks=t_vals[1]:2.0:t_vals[end]
)
for i in 1:N
    lines!(ax2, t_vals, getindex.(v_vals, i), label=L"\dot{q}_{%$i}")  # Velocities
    #scatter!(ax2, T, ywant[N+i, :], color=:black, marker=:x, markersize=8, label=L"$x$")  # Desired 
end
axislegend(ax2, unique=true, orientation=:horizontal)
#
base_vals = zeros(o, 2N, length(t_vals))
# 3. Base functions plot using pi_paper
ax3 = Axis(f[3, 1],
    title="Base Functions",
    xlabel="Time",
    ylabel=L"$\pi_{i}$",
    ylabelrotation=0,  # Rotates y-label 90 degrees
    ylabelsize=18,
    xticks=t_vals[1]:2.0:t_vals[end]
)
#= for i in 1:o
    base_vals[i, :, :] .= hcat([pi_paper(A, B, C, t, T, i) for t in t_vals]...)
    for j in 1:2N# Using pi_paper to compute the base functions
        color_index = ((i - 1) * 2N + j) / (o * 2N)
        lines!(ax3, t_vals, base_vals[i, j, :], color=plotcolors[color_index])
    end
end =#
#
# 4. Control Action plot (input uₙ)
ax4 = Axis(f[4, 1],
    title="Control Action",
    xlabel="Time",
    ylabel=L"$u(t)$",
    ylabelrotation=0,  # Rotates y-label 90 degrees
    ylabelsize=18,
    xticks=t_vals[1]:2.0:t_vals[end]
)
#lines!(ax4, t_vals, u_vals, color=:black)

# Display the figure
display(f)