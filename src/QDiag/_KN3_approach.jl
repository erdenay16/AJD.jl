include("_utils.jl")

function _optimize(
    C_0,
    C,
    weights,
    approach,
    tolerance,
    maximum_iteration,
    random_number_generator,
)
    if approach == "KN3"
        result, iteration_errors = _apply_KN3(
            C_0,
            C,
            weights,
            tolerance,
            maximum_iteration,
            random_number_generator,
        )
    else
        # N5 will be implemented
    end
    return result, iteration_errors
end

function _initialize_KN3(C_0, C, weights, random_number_generator)
    N = size(C_0, 1)
    P = _find_P(C_0)
    W = randn(random_number_generator, size(C_0))
    K = size(C, 3)

    for k = 1:K
        C[:, :, k] = P' * C[:, :, k] * P
    end

    M = zeros(N, N)
    weights = weights ./ norm(weights, 1)
    return P, W, C, M, weights
end

function _M_loop_KN3(M, C, W, weights)
    K = size(C, 3)
    for k = 1:K
        M_1 = C[:, :, k] * W'
        M_2 = C[:, :, k]' * W'
        M = M + weights[k] * (M_1 * M_1' + M_2 * M_2')
    end
    return M
end

function _calculate_error(W, C_0, C, weights)
    cumulative_error = 0
    N = size(W, 1)
    d = diag(W * C_0 * W')
    W = W ./ sqrt.(d')
    D = permutedims(C, (2, 1, 3)) .* W
    D = permutedims(D, (2, 1, 3)) .* W
    D = D.^2

    for i in 1:length(weights)
        cumulative_error += sum(D[:, :, i] - Diagonal(D[:,:,i]))
    end
    return cumulative_error / (N^2 - N)
end

function _apply_KN3(C_0, C, weights, tolerance, maximum_iteration, random_number_generator)
    P, W, C, M, weights = _initialize_KN3(C_0, C, weights, random_number_generator)
    M = _M_loop_KN3(M, C, W, weights)
    N = size(C_0, 1)
    K = size(C, 3)
    iterations = 1
    delta = 0.0
    iteration_errors = Float64[]

    while iterations <= maximum_iteration
        for i = 1:N
            w_i = W[i, :]
            for k = 1:K
                m_1 = C[:, :, k] * w_i
                m_2 = C[:, :, k]' * w_i
                M = M - weights[k] * (m_1 * m_1' + m_2 * m_2')
            end
            eigen_values_M, eigen_vectors_M = eigen(M)
            w_new = eigen_vectors_M[:, argmin(eigen_values_M)]
            for k = 1:K
                m_1 = C[:, :, k] * w_new
                m_2 = C[:, :, k]' * w_new
                M = M + weights[k] * (m_1 * m_1' + m_2 * m_2')
            end

            delta = norm(W[i, :] - w_new)
            W[i, :] = w_new
        end

        push!(iteration_errors, _calculate_error(W, C_0, C, weights))
        if delta <= tolerance
            break
        end

        iterations += 1
    end
    return W * P', iteration_errors
end
