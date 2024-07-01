using LinearAlgebra: Diagonal, diag, diagm
using Random: randn
using Statistics: std


function genFullyDiagMs(K::Int, size::Int)

    C0 = Diagonal(ones(size))
    Cdiag = [Diagonal(randn(size)) for _ in 1:K]
    Cs = pushfirst!(Cdiag, C0)

    A = randn(size, size)

    Cx = Array{Float64,3}(undef, size, size, K + 1)
    for i in 1:K+1
        Cx[:, :, i] = A * Cs[i] * A'
    end
    return Cx, A
end


function genApproxDiagMs(K::Int, size::Int)
    C0 = Diagonal(ones(size))

    Cdiag = [Diagonal(randn(size)) for _ in 1:K]
    Cs = pushfirst!(Cdiag, C0)

    A = randn(size, size)
    Ax = [A]

    Cx = Array{Float64,3}(undef, size, size, K + 1)
 
    for i in 1:K+1
        Ak = A .+ randn(size, size) .* (10^(-2))
        Cx[:, :, i] = Ak * Cs[i] * Ak'
        push!(Ax, Ak)
    end 
    return Cx, Ax  
end

function check_diag(Cs)
    tol = 1e-9 
    _, _, K = size(Cs)
    for k in 1:K
        off_diag_sum = sum(abs, Cs[:, :, k]) - sum(abs, diagm(diag(Cs[:, :, k])))
        return off_diag_sum < tol
    end
end