using LinearAlgebra: I, diag, diagm, norm, tr, opnorm
using Base: frexp


"""
ffdiag(Cs, runs, tol) -> AbstractArray{<:Real}, Matrix{<:Real}, Array{<:Real}

This is an Implementation of the Algorithm introduced in: 
Ziehe, Andreas; Laskov, Pavel; Nolte, Guido; Müller Klaus-Robert. (2004).
A Fast Algorithm for Joint Diagonalization with Non-orthogonal
Transformations and its Application to
Blind Source Separation.
Journal of Machine Learning Research 5 (2004) 777–800.

The function returns the diagonalized set of matrices, the diagonalization matrix and an array of diagonlization errors per iteration.

# Arguments
- Cs::AbstractArray{<:Real}: This is a set of matrices to be diagonalized.
- runs::Int: The maximum number of iterations. The default is max 100 iterations.
- tol::Float64: The tolerance for the error. The default is 1e-9.
"""

function ffdiag(
    C0::AbstractArray{<:Real},
    runs::Int=100,
    tol::Float64=1e-9
)
    @assert length(size(C0)) == 3 "C must be a tensor with dimensions K x M x N"
    @assert tol >= 0 "Tolerance must be a nonnegative real number."
    @assert runs > 0 "Maximum iteration must be positive."

    dim1, dim2, K = size(C0) # m,n,K

    @assert (dim1 == dim2) "Error: Matrices not square."
    @assert K > 1 "Error: Input is only one matrix not a set of matrices."

    # Assert symmetry of matrices
    for k in 1:K
        @assert isapprox(C0[:, :, k], C0[:, :, k]') "Error: Matrix Cs[$k] is not symmetrical."
    end

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


# W is calculated as described by equation 17 in the article. 

function getW(Cs)

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



# Returns the magnitude of the off-diagonal elements. When f apporaches zero the matrix approaches diagonalization

function off(V, C)
    F = V * C * V'
    f = tr(F' * F) - tr(F .* F)
    return f 
end


# Returns the sum of the magnitude of the off-diagonal elements for multiple matrices 

function get_off(V, Cs)
    _, _, K = size(Cs)
    f = 0

    for k in 1:K
        f = f + off(V, Cs[:, :, k])
    end
    return f
end


# Calculates the norm of the elements that are not on the diagonal and sums it over all matirces of C
function cost_off(Cs, V)
    n, m, K = size(Cs)

    cost = 0
    for k in 1:K
        Ck = (V * Cs[:, :, k] * V')
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

# norms the the matrix 
function normit(V)
    N, _ = size(V)
    for n in 1:N
        nn = norm(V[:, n])
        if nn >= 1e-9
            V[:, n] = V[:, n] ./ nn
        end
    end
    return V
end