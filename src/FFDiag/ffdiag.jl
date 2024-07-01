
using LinearAlgebra: I, diag, diagm, norm, tr, opnorm
using Base: frexp


"""
ffdiag(Cs) -> AbstractArray{<:Real}

This is an Implementation of the Algorithm introduced in: 
Ziehe, Andreas; Laskov, Pavel; Nolte, Guido; Müller Klaus-Robert. (2004).
A Fast Algorithm for Joint Diagonalization with Non-orthogonal
Transformations and its Application to
Blind Source Separation.
Journal of Machine Learning Research 5 (2004) 777–800.

# Arguments
- Cs::AbstractArray{<:Real}: This is a set of matrices to be diagonalized.
"""

function ffdiag(
    C0::AbstractArray{<:Real},
    runs::Real=100,
    tol::Float64=1e-9 # = eps
)

    @assert length(size(C0)) == 3 "C must be a tensor with dimensions K x M x N"
    @assert tol >= 0 "Tolerance must be a nonnegative real number."
    @assert runs > 0 "Maximum iteration must be positive."

    dim1, dim2, K = size(C0) # m,n,K

    @assert (dim1 == dim2) "Error: Matrices not square"

    # Init
    V = I
    Cs = copy(C0)

    df = 1
    run = 1
    errs = zeros(runs + 1)
    fs = zeros(runs + 1)

    # doesnt do anything right now because V = I 
    for k in 1:K
        V * Cs[k] * V'
    end

    while run < runs && df > tol

        W = getW(Cs)

        _, e = frexp(opnorm(W, Inf))

        s = max(0.0, e - 1)
        W = W / (2^s)

        # Compute V
        V = (I + W) * V

        # Renormalization
        V = diagm(1 ./ sqrt.(diag(V * V'))) * V

        for k in 1:K
            # Cs[:, :, k] = ((I + W) * C0[:, :, k] * (I + W)')
            Cs[:, :, k] = (V * C0[:, :, k] * V')
            #errs[run+1] += norm(Cs[:, :, k] - diagm(diag(Cs[:, :, k]))) / K
        end

        fs[run] = get_off(V, C0)
        errs[run] = cost_off(C0, normit(V')')

        if run > 2
            df = abs(fs[run-1] - fs[run])
        end

        run += 1
    end

    return Cs, V, errs[1:run-1]
end



function getW(Cs::AbstractArray{<:Real})

    dim1, dim2, K = size(Cs)

    W = zeros(dim1, dim2)
    Ds = [zeros(dim1) for _ in 1:K]
    Es = copy(Cs)

    for i in 1:K
        Ds[i] = (diag(Cs[:, :, i]))
        # make the diagonal of Es zero
        for j in 1:dim1
            Es[j, j, i] = 0
        end
    end

    # calculate W (17) 
    z = zeros(dim1, dim2)
    y = zeros(dim1, dim2)

    for k in 1:K
        z += (Ds[k] * Ds[k]')
        y += 0.5 .* Ds[k]' .* (Es[:, :, k] + Es[:, :, k]')
    end

    # Check if we need to the for loops are programmed correctly 
    for i in 1:dim1-1
        for j in (i+1):dim2
            if i != j
                W[i, j] = (z[i, j] * y[j, i] - z[i, i] * y[i, j]) / (z[j, j] * z[i, i] - z[i, j]^2)
                W[j, i] = (z[i, j] * y[i, j] - z[j, j] * y[j, i]) / (z[j, j] * z[i, i] - z[i, j]^2)
            end
        end
    end
    return W
end


function off(V, C)
    F = V * C * V'
    f = tr(F' * F) - tr(F .* F)
    return f 
end

function get_off(V, C)
    _, _, K = size(C)
    f = 0

    for k in 1:K
        f = f + off(V, C[:, :, k])
    end
    return f
end

function cost_off(C, V)
    n, m, K = size(C)

    cost = 0
    for k in 1:K
        Ck = (V * C[:, :, k] * V')
        for i in 1:n
            for j in 1:m
                if i == j
                    Ck[i, j] = 0
                end
            end
        end
        cost += norm(Ck)^2
    end
    return cost

end

function normit(V)
    N, M = size(V)
    V_res = copy(V)
    for n in 1:N
        nn = norm(V[:, n])
        if nn >= 1e-9
            V[:, n] = V[:, n] ./ nn
        end
    end
    return V_res
end
