using LinearAlgebra: I, diag, diagm, norm, tr, opnorm
using Base: frexp


"""
    ffdiag(C0, runs=100, tol=1e-9) -> AbstractArray{<:Real}, Matrix{<:Real}, Array{<:Real}

Compute the transformation matrix that diagonalizes a set of symmetric matrices

This is an Implementation of the Algorithm introduced in: 
Ziehe, Andreas; Laskov, Pavel; Nolte, Guido; Müller Klaus-Robert. (2004).
A Fast Algorithm for Joint Diagonalization with Non-orthogonal
Transformations and its Application to
Blind Source Separation.
Journal of Machine Learning Research 5 (2004) 777–800.

The function returns the diagonalized set of matrices, the diagonalization matrix and an array of diagonlization errors per iteration.

# Arguments
- C0::AbstractArray{<:Real}: This is a set of matrices to be diagonalized.
- runs::Int: The maximum number of iterations. The default is max 100 iterations.
- tol::Float64: The tolerance for the error. The default is 1e-9.
"""
function ffdiag(
    C0::AbstractArray{<:Real},
    runs::Int=100,
    tol::Float64=1e-9
)
    length(size(C0)) == 3 || throw(DimensionMismatch("C must be a tensor with dimensions K x M x N"))
    tol >= 0 || throw(ArgumentError("Tolerance must be a nonnegative real number."))
    runs > 0 || throw(ArgumentError("Maximum iteration must be positive."))

    dim1, dim2, K = size(C0) # m,n,K

    (dim1 == dim2) || throw(DimensionMismatch("Error: Matrices not square."))
    K > 1 || throw(ArgumentError("Error: Input is only one matrix not a set of matrices."))

    # Assert symmetry of matrices
    for k in 1:K
        isapprox(C0[:,:,k], C0[:,:,k]') || @warn "Error: Matrix Cs[$k] is not symmetrical."
    end

    # Init
    V = I
    Cs = copy(C0)

    df = 1
    run = 1
    errs = zeros(runs + 1)
    fs = zeros(runs + 1)

    while run < runs && df > tol

        W = getW(Cs)

        _, e = frexp(opnorm(W, Inf))

        s = max(0.0, e - 1)
        W = W / (2^s)

        # Compute V
        V = (I + W) * V

        # Renormalization
        V = diagm(1 ./ sqrt.(diag(V * V'))) * V

        # Update Cs
        [Cs[:, :, k] = (V * C0[:, :, k] * V') for k in 1:K]

        # Compute off-diagonal elements 
        fs[run]   = get_off(V, C0)
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
    
    # Ds = the diagonal of Cs
    Ds = [diag(C) for C in eachslice(Cs, dims=3)]
    
    # Es = Cs with the diagonal set to zero
    Es = copy(Cs); [Es[:, :, i] -= diagm(Ds[i]) for i in 1:K]
    
    # calculate W (17) 
    z = zeros(dim1, dim2); 
    y = zeros(dim1, dim2)
    
    for k in 1:K
        z += (Ds[k] * Ds[k]')
        y += 0.5 .* Ds[k]' .* (Es[:, :, k] + Es[:, :, k]')
    end
    
    W = zeros(dim1, dim2)

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
    f = sum(off(V, C) for C in eachslice(Cs, dims=3))
    return f
end


# Calculates the norm of the elements that are not on the diagonal and sums it over all matirces of C
function cost_off(Cs, V)
    K = size(Cs, 3)

    cost = 0
    for C in eachslice(Cs, dims=3)
        Ck = (V * C * V')
        # The diagonal is set to zero
        Ck -= diagm(diag(Ck))
        cost += norm(Ck)^2
    end
    return cost

end

# norms the the matrix 
function normit(V)
    for col in eachcol(V)
        nn = norm(col)
        if nn >= 1e-9
            col .= col / nn
        end
    end
    return V
end