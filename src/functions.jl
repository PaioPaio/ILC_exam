# A lot of repeated code but whatever
"""
Underactuated Mass-Spring-Damper with N carts. The actuation appears as a Force on the last cart.
"""
function MSD_matrices(N, m, κ, σ)
    A = zeros(2N, 2N)
    B = zeros(2N, 1)

    A[1:N, N+1:2N] = I(N)  # dx/dt = v
    A[N+1, 1:2] = [-2κ / m, κ / m]
    A[N+1, N+1:N+2] = [-2σ / m, σ / m]
    for i in N+2:2N-1
        A[i, i-N-1:i-N+1] = [κ / m, -2κ / m, κ / m]
        A[i, i-1:i+1] = [σ / m, -2σ / m, σ / m]
    end
    A[2N, N-1:N] = [κ / m, -κ / m]
    A[2N, 2N-1:2N] = [σ / m, -σ / m]

    B[2N] = 1 / m

    C = float(I(2N))

    return A, B, C
end
"""
Basketball in the wind system. The gravity constant is considered unknown. The basketball player can act only whenever the ball is touching the ground
"""
function basket_matrices(m)
    A = zeros(6, 6)
    B = zeros(6, 3)
    C = zeros(3, 6)

    A[1:3, 4:6] = I(3)
    B[4:6, :] = I(3) ./ m

    C[:, 1:3] = I(3)

    return A, B, C
end
"""
Convenience function to initialize containers of the solution
"""
function initialize_sol(Ns, No, l, o, Niter)
    x = zeros(div(Ns, 2), o, Niter)
    v = zeros(div(Ns, 2), o, Niter)
    u = zeros(l, o, Niter)
    y = zeros(No, o, Niter)
    αrray = zeros(No, o, Niter)
    solarray = Vector{Any}(undef, Niter)

    return x, v, u, y, αrray, solarray
end


"""
Compute H matrix by integrating Grahamian Matrices
"""
function compute_H(A, B, C, T, fun)
    Twith0 = [0.0; T]
    No, Ns = size(C)
    o = size(T, 1)
    H = zeros(No * o, No * o)
    for i in 1:o, j in 1:o
        if i >= j
            # f(t,p) where t is the integration variable and p is tuple of parameters (we don't need that since we are directly giving the function A,B,C)
            f(t, p) = fun(A, B, C, Twith0[i+1] - t) * fun(A, B, C, Twith0[j+1] - t)'
            H[(i-1)*No+1:i*No, (j-1)*No+1:j*No] .= solve(IntegralProblem(f, (Twith0[j], Twith0[j+1])), HCubatureJL()).u
        else
            H[(i-1)*No+1:i*No, (j-1)*No+1:j*No] .= zeros(No, No)
        end
    end
    return H
end

compute_H(A, B, C, T) = compute_H(A, B, C, T, CAtB)

function compute_H_constr(A, B, C, T, dT, fun)
    Twith0 = [0.0; T]
    No, Ns = size(C)
    o = size(T, 1)
    H = zeros(No * o, No * o)
    for i in 1:o, j in 1:o
        if i >= j
            # f(t,p) where t is the integration variable and p is tuple of parameters (we don't need that since we are directly giving the function A,B,C)
            f(t, p) = fun(A, B, C, Twith0[i+1] - t) * fun(A, B, C, Twith0[j+1] - t)'
            #Main.@infiltrate
            H[(i-1)*No+1:i*No, (j-1)*No+1:j*No] .= solve(IntegralProblem(f, (max(0.0, Twith0[j]), min(Twith0[j+1], Twith0[j] + dT))), HCubatureJL()).u
        else
            H[(i-1)*No+1:i*No, (j-1)*No+1:j*No] .= zeros(No, No)
        end
    end
    return H
end

compute_H_constr(A, B, C, T, dT) = compute_H_constr(A, B, C, T, dT, CAtB)

"""
Computes the basis functions used in the paper  given the state matrices
"""
function pi_paper(A, B, C, t, T, i)
    No, Ns = size(C)
    l = size(B, 2)
    minT = [0.0; T[1:end-1]]
    if t < T[i] && t >= minT[i]
        return CAtB(A, B, C, T[i] - t)
    else
        return zeros(No, l)
    end
end

function CAtB(A, B, C, t)
    return C * exp(A * t) * B
end

"""
Compute input functions given a basis and an alpha vector
"""
function utotn(A, B, C, α, f, t, T)
    l = size(B, 2)
    minT = [0.0; T[1:end-1]]
    return sum(
        if t < T[i] && t >= minT[i]
            f(A, B, C, T[i] - t)' * α[:, i]
        else
            zeros(l)
        end
        for i in eachindex(T)
    )
end
# Method with standard basis
utotn(A, B, C, α, t, T) = utotn(A, B, C, α, CAtB, t, T)

"""
Constrained Version of utotn
"""
function utotconstr(A, B, C, α, f, t, T, dT)
    l = size(B, 2)
    Twith0 = [0.0; T]
    return sum(
        if t < Twith0[i] + dT && t >= Twith0[i]
            Main.@infiltrate
            f(A, B, C, T[i] - t)' * α[:, i]
        else
            zeros(l)
        end
        for i in eachindex(T)[1:end-1]
    )
end
# Method with standard basis
utotconstr(A, B, C, α, t, T, dT) = utotconstr(A, B, C, α, CAtB, t, T, dT)

"""
Linear dynamics, either with u (f in the body function, while u is the state just because the ODE library we're using has this convention) as a parameter or directly coded inside the function.
"""
function LinearDynamics!(du, u, p, t)
    A, B, C, α, T, f = p

    du .= A * u + B * f(A, B, C, α, t, T)
end

function LinearDynamics_intrinsic!(du, u, p, t)
    x, v, input = u[1:N], u[N+1:2N], u[end]
    A, B, C, α, T = p

    du[1:2N] = A * u[1:2N] + B * input[1]
    du[end] = utotn(A, B, C, α, t, T) - input[1]
end

"""
Linear Dynamics with constant unknown offset on inputs.
"""
function LinearDynamicsDisturbed!(du, u, p, t)

    A, B, C, α, T, dT, f, dist = p

    du .= A * u + B * (f(A, B, C, α, t, T, dT) + dist)

end