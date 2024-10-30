using OrdinaryDiffEq, LinearAlgebra, Integrals, CairoMakie, ColorSchemes

includet("../src/functions.jl")

# System parameters
N = 7
T = [8.0, 18.0, 20.0]
o = length(T)
m, κ, σ = 1.0, 2.0, 1.0
pp = (m, κ, σ)
Niter = 30

# Define system matrices

A, B, C = MSD_matrices(N, m, κ, σ)
##
# Desired output
x₁ = cumsum(ones(N)) ./ N
y₁ = vcat(x₁, zeros(N))
ywant = hcat(y₁, zeros(2N), zeros(2N))

# Initialize storage for results
xⁱ = zeros(N, o, Niter)
vⁱ = zeros(N, o, Niter)
uⁱ = zeros(1, o, Niter)
yⁱ = zeros(2N, o, Niter)
αrray = zeros(2N, o, Niter)
solarray = Vector{Any}(undef, Niter)

##  
# Compute H matrix

H = compute_H(A, B, C, T)               # Lower block triangular as axpected
L = inv(H'H + I(2N * o) * 1e-2) * H'    # Levenberg-Marquaad (or something along that way) with ϵ = 1e-2

# ILC loop
u0 = vcat(zeros(N), zeros(N), 0.0)
tspan = (0.0, T[end])
for i in 1:Niter
    # Set up and solve ODE
    p = (A, B, C, αrray[:, :, i], T, utotn)
    prob = ODEProblem(LinearDynamics!, u0, tspan, p)
    sol = solve(prob, Tsit5(), abstol=1e-9)
    solarray[i] = sol

    # Extract and store results
    for (j, t) in enumerate(T)
        xⁱ[:, j, i] = sol(t)[1:N]
        vⁱ[:, j, i] = sol(t)[N+1:2N]
    end

    # Update α for next iteration
    yⁱ[:, :, i] = vcat(xⁱ[:, :, i], vⁱ[:, :, i])
    if i < Niter
        αrray[:, :, i+1] = αrray[:, :, i] + reshape(L * (ywant[:] - yⁱ[:, :, i][:]), 2N, o)
    end
end

##
# Time vector from the solution (full time points)
t_vals = 0.0:0.05:T[end]
x_vals = [solarray[end](t)[1:N] for t in t_vals]  # Position data
v_vals = [solarray[end](t)[N+1:2N] for t in t_vals]  # Velocity data 
u_vals = [utotn(A, B, C, αrray[:, :, end], t, T) for t in t_vals]
errors = norm.([(ywant[:] - yⁱ[:, :, i][:]) for i in 1:Niter])

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
    scatter!(ax2, T, ywant[N+i, :], color=:black, marker=:x, markersize=8, label=L"$x$")  # Desired 
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
for i in 1:o
    base_vals[i, :, :] .= hcat([pi_paper(A, B, C, t, T, i) for t in t_vals]...)
    for j in 1:2N# Using pi_paper to compute the base functions
        color_index = ((i - 1) * 2N + j) / (o * 2N)
        lines!(ax3, t_vals, base_vals[i, j, :], color=plotcolors[color_index])
    end
end
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
lines!(ax4, t_vals, u_vals, color=:black)

display(f)

## Save figure 
imgpath = "./report_slides/typst/images/"
figname = "MSD_$(N)carts_$(Niter)iters.svg"
save(imgpath * figname, f)