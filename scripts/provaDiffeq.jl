using OrdinaryDiffEq, LinearAlgebra, Integrals, GLMakie

# System parameters
const N = 5
const o = 3
const T = [8.0, 18.0, 20.0]
const m, κ, σ = 1.0, 2.0, 1.0
const Niter = 10

# Define system matrices
function system_matrices(N, m, κ, σ)
    A = zeros(2N, 2N)
    B = zeros(2N, 1)
    
    A[1:N, N+1:2N] = I(N)  # dx/dt = v
    A[N+1, 1:2] = [-2κ/m, κ/m]
    for i in N+2:2N-1
        A[i, i-N-1:i-N+1] = [κ/m, -2κ/m, κ/m]
        A[i, i-1:i+1] = [σ/m, -2σ/m, σ/m]
    end
    A[2N, N-1:N] = [κ/m, -κ/m]
    A[2N, 2N-1:2N] = [σ/m, -σ/m]
    
    B[N+1] = 1/m
    B[2N] = 1/m
    
    return A, B
end

A, B = system_matrices(N, m, κ, σ)
C = float(I(2N))  # y = [x; v], so C is identity

function CAtB(A, B, C, t)
    return C * exp(A * t) * B
end

function utotn(A, B, C, α, t, T)
    minT = [0.0; T[1:end-1]]
    return sum(
        if t <= T[i] && t > minT[i]
            dot(CAtB(A, B, C, T[i]-t), α[:,i])
        else
            0.0
        end
        for i in eachindex(T)
    )
end

function upi(A,B,C,t,T,i)
    minT = [0.0; T[1:end-1]]
    if t <= T[i] && t > minT[i]
        return CAtB(A, B, C, T[i]-t)'
    else
        return 0.0
    end
end

function mass_spring_damper!(du, u, p, t)
    x, v, input = u[1:N], u[N+1:2N], u[end]
    A, B, C, α, T = p

    du[1:2N] = A * u[1:2N] + B * input[1]
    du[end] = utotn(A, B, C, α, t, T) - input[1]
end

# Desired output
x₁ = cumsum(ones(N) ./ N)
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
function compute_H(A, B, C, T)
    Twith0 = [0.0; T]
    H = zeros(2N*o, 2N*o)
    for i in 1:o, j in 1:o
        if i >= j
            f(t,p) = C * exp(A * (Twith0[i+1] - t)) * B * B' * exp(A' * (Twith0[j+1] - t)) * C'
            result = solve(IntegralProblem(f,(Twith0[j],Twith0[j+1])),HCubatureJL()).u
            H[(i-1)*2N+1:i*2N, (j-1)*2N+1:j*2N] .= result
        else
            H[(i-1)*2N+1:i*2N, (j-1)*2N+1:j*2N] .= zeros(2N,2N)
        end
    end
    return H
end

H = compute_H(A, B, C, T)
L = inv(H'H + I(2N*o)*1e-2) * H'

## ILC loop
for i in 1:Niter
    # Set up and solve ODE
    u0 = vcat(zeros(N), zeros(N), 0.0)
    tspan = (0.0, T[end])
    p = (A, B, C, αrray[:,:,i], T)
    prob = ODEProblem(mass_spring_damper!, u0, tspan, p)
    sol = solve(prob, Tsit5())
    solarray[i] = sol

    # Extract and store results
    for (j, t) in enumerate(T)
        xⁱ[:, j, i] = sol(t)[1:N]
        vⁱ[:, j, i] = sol(t)[N+1:2N]
        uⁱ[:, j, i] = sol(t)[end:end]
    end

    # Update α for next iteration
    if i < Niter
        y_current = vcat(xⁱ[:,:,i], vⁱ[:,:,i])
        αrray[:,:,i+1] = αrray[:,:,i] + reshape(L * (ywant[:] - y_current[:]), 2N, o)
    end
end

##
# Time vector from the solution (full time points)
t_vals = solarray[end].t  
x_vals = [solarray[end](t)[1:N] for t in t_vals]  # Position data
v_vals = [solarray[end](t)[N+1:2N] for t in t_vals]  # Velocity data
u_vals = [solarray[end](t)[end] for t in t_vals]  # Control input over time

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

# 3. Base functions plot using upi
#= ax3 = Axis(f[3, 1], title="Base Functions π", xlabel="Time", ylabel="Base Functions")
for i in 1:o
    base_vals = [upi(A, B, C, t, T, i) for t in t_vals]  # Using upi to compute the base functions
    lines!(ax3, t_vals, base_vals, color=RGB(0.1*i, 0.8-i/10, 0.9-i/10))
end =#

# 4. Control Action plot (input uₙ)
ax4 = Axis(f[4, 1], title="Control Action uₙ", xlabel="Time", ylabel="Control Action")
lines!(ax4, t_vals, u_vals, color=:black)

# Display the figure
display(f)


##

# Set up the animation
scene = Scene(resolution = (800, 400))

# Set up the masses (represented as circles) and springs (lines between them)
mass_radius = 0.5  # The size of the carts
wall_position = -1.0  # The position of the wall
t_vals = solarray[end].t  # Continuous time points from the last iteration
x_vals = [solarray[end](t)[1:N] for t in t_vals]  # Position data over time

# Initialize the mass graphics
masses = [scatter!([x_vals[1][i]], [0.0], markersize=mass_radius*100, color=:blue) for i in 1:N]

# Initialize the spring graphics (connecting lines between masses)
springs = []
for i in 1:N-1
    push!(springs, lines!([x_vals[1][i], x_vals[1][i+1]], [0.0, 0.0], linewidth=2, color=:black))
end
# Add a line for the wall-to-first-mass spring
push!(springs, lines!([wall_position, x_vals[1][1]], [0.0, 0.0], linewidth=2, color=:black))

# Function to update the positions of masses and springs in the animation
function update_scene(timestep)
    for i in 1:N
        # Clear previous positions
        scatter!(masses[i], [x_vals[timestep][i]], [0.0], markersize=mass_radius*100, color=:blue)
    end
    # Update the springs (lines between carts)
    for i in 1:N-1
        lines!(springs[i], [x_vals[timestep][i], x_vals[timestep][i+1]], [0.0, 0.0], linewidth=2, color=:black)
    end
    # Update the wall-to-first-mass spring
    lines!(springs[N], [wall_position, x_vals[timestep][1]], [0.0, 0.0], linewidth=2, color=:black)
end

# Record the animation and save it
record(scene, "mass_spring_animation.mp4", 1:length(t_vals); framerate=30) do timestep
    update_scene(timestep)
end
