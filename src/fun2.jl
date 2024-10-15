function system_matrices(N, m, κ, σ)
    A = zeros(2N, 2N)
    B = zeros(2N, 1)
    
    A[1:N, N+1:2N] = I(N)  # dx/dt = v
    A[N+1, 1:2] = [-2κ/m, κ/m]
    A[N+1, N+1:N+2] = [-2σ/m, σ/m]
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
        return CAtB(A, B, C, T[i]-t)    
    else
        return zeros(2N)
    end
end

function mass_spring_damper2!(du, u, p, t)
    A, B, C, α, T, f = p

    du[1:2N] = A * u[1:2N] + B * f(A, B, C, α, t, T)
end

function mass_spring_damper!(du, u, p, t)
    x, v, input = u[1:N], u[N+1:2N], u[end]
    A, B, C, α, T = p

    du[1:2N] = A * u[1:2N] + B * input[1]
    du[end] = utotn(A, B, C, α, t, T) - input[1]
end

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