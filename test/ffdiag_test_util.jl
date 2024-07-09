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


# function genTimeCorrMs(X, lags, iscov=true, symmetrize=false)
#     K = length(lags)

#     T, N = size(X)

#     if T < N
#         @info "X may be oriented the wrong way" 
#     end

#     for i in 1:k
#         if lags(k) > 0
#             x1 = x[1:end-lags[k], :]
#             x2 = x[1+lags[k]:end, :]
#         else
#             x1 = [1-lags[k]:end, :]
#             x2 = [1:end+lags[k], :]
#         end

#         T = size(x1, 1)

#         if iscov == true
#             x1 = x1 .- repeat(sum(x1, dims=1) / T, T, 1)
#             x2 = x2 .- repeat(sum(x2, dims=1) / T, T, 1)
#             C[:, :, k] = x1' * x2 / (T-1)
#         else
#         end
    
#     end
# end
