# N5 Approach

using LinearAlgebra
using Random

function _initialize_N5(C_0, C, weights, random_number_generator)
    N = size(C_0, 1)
    W = randn(random_number_generator, N, N)
    P = cholesky(C_0).U
    for k = 1:size(C, 3)
        C[:, :, k] = P' * C[:, :, k] * P
    end
    weights = weights / sum(weights)
    return P, W, C, weights
end

function _optimize_N5(C, W, weights, tolerance, max_iter)
    N = size(C, 1)
    K = size(C, 3)
    for iteration in 1:max_iter
        M = zeros(N, N)
        for k = 1:K
            W_C_k = W' * C[:, :, k] * W
            M += weights[k] * (W_C_k * W_C_k' + W_C_k' * W_C_k)
        end
        
        eigen_values, eigen_vectors = eigen(M)
        W_new = eigen_vectors[:, argmin(eigen_values)]
        
        error = norm(W - W_new)
        W = W_new
        
        if error < tolerance
            break
        end
    end
    return W
end

function QDiag_N5(C_0::AbstractArray{<:Real}, C::AbstractArray{<:Real}, weights::AbstractArray{<:Real}, tolerance::Real, max_iter::Int, rng::Random.AbstractRNG = Random.default_rng())
    @assert isposdef(C_0) "C_0 must be positive-definite."
    @assert size(C_0, 1) == size(C_0, 2) == size(C, 1) == size(C, 2) "C_0 and C must have compatible dimensions."
    @assert length(size(weights)) == 1 "Weights must be a vector."
    @assert size(weights, 1) == size(C, 3) "Weights vector length must match the third dimension of C."
    
    P, W, C, weights = _initialize_N5(C_0, C, weights, rng)
    W = _optimize_N5(C, W, weights, tolerance, max_iter)
    return P * W
end