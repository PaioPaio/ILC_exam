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

    function expnumsymb(A::AbstractMatrix{T},t::Num) where T 
        evals, evec = eigen(A)
        return evec * diagm(exp.(evals.*t)) * evec^-1
    end


    expsym(A::AbstractMatrix{Num},t::T) where T = expmapproxint(A,t)
    expsym(A::AbstractMatrix{T},t::Num) where T = expnumsymb(A,t)
    expsym(A::AbstractMatrix{Num},t::Num) = expmapproxint(A,t)
    expsym(A::AbstractMatrix{T},t::T) where T = exp(A*t) 

    ## Function Basis
    function CAtB(A,B,C,t)
        return C * expsym(A,t) * B
    end
    CAtB(A::AbstractMatrix{T},B::AbstractMatrix{T},C::AbstractMatrix{T},t::T) where T = C*exp(A*t)*B

    function pif(A,B,C,T,i,t)
        minT = i==1 ? 0 : T[i-1]
        maxT = T[i]
        return ifelse(!((t>maxT)&(t<minT)), CAtB(A,B,C,maxT-t)', zeros(size(B,2),size(C,1)))
        #return ifelse(!((t>T[i])&(t<T[i-1])), B' * (exponential(A'*(T[i]-t))) * C', zeros(size(B,2),size(C,1)))
    end

    function piall(A,B,C,T,t)

        minT = [0.0 T[1:end-1]...]
        return hcat(ifelse(!((t>T[i])&(t<minT[i])), CAtB(A,B,C,T[i]-t), 0.0) for i in eachindex(T))

    end

    function uiter(A,B,C,T,t,α)
        minT= [0.0 T[1:end-1]...]       
        return sum(ifelse(!((t>T[i])&(t<minT[i])), dot(CAtB(A,B,C,T[i]-t),α[:,i]), 0.0) for i in eachindex(T))
    end

    function cleanupimaginary(arraytoclean::AbstractVecOrMat{Num};tol=1e-8)
        indexes = map(x->Symbolics.coeff(arraytoclean))
    end

    function cleanupimaginary(arraytoclean::AbstractVecOrMat{T};tol=1e-8) where T
    end