using LinearAlgebra

function fast_frobenius(Cs)
    """
    Ziehe (2004):
    A Fast Algorithm for Joint Diagonalization with Non-orthogonal
    Transformations and its Application to
    Blind Source Separation
    
    Input   Cs  Matrices to be diagonalized
    Output  V diagonalizer that minimizes off diagonal terms of Cs
            errs average error
    """
    sweeps=100
    theta = 0.5
    K, m, N = size(Cs)
    # @assert m == N
    W = zeros(m, m)
    _, V = eigen(Cs[:,:,1])  # first guess for diagonalizer

    z = zeros(m, m)
    y = zeros(m, m)
    
    Ds = zeros(K, m)   # diagonal terms of Cs
    Es = similar(Cs)  # offdiagonal terms of Cs
    
    errs = zeros(sweeps+1)
    
    for C in eachslice(Cs, dims=3)
        # calculate average initial error
        C = V * C * V'
        errs[1] += norm(C - Diagonal(diag(C)))/K
    end

    for s in 1:3
        for k in 1:K
            # set Ds and Es
            Ds[k, :] = diagm(Cs[:, :, k])
            Es[k, :, :] = Cs[:, :, k]
            fill!(Es[k, :, :], 0)
        end
        z = zeros(m, m)
        y = zeros(m, m)
        
        # compute W from Cs according to equation 17 in article
        for i in 1:m
            for j in 1:m
                for k in 1:K
                    z[i, j] += Ds[k, i] * Ds[k, j]
                    y[i, j] += 0.5 * Ds[k, j] * (Es[k, i, j] + Es[k, j, i])
                end
            end
        end

        for i in 1:m-1
            for j in i+1:m
                W[i, j] = (z[j, i] * y[j, i] - z[i, i] * y[i, j]) / (z[j, j] * z[i, i] - z[i, j] * z[i, j])
                W[j, i] = (z[i, j] * y[i, j] - z[j, j] * y[j, i]) / (z[j, j] * z[i, i] - z[i, j] * z[i, j])
            end
        end

        # make sure W satisfies frobenius norm < theta
        while norm(W, "fro") > theta
            W = W * (theta / norm(W, "fro"))
        end

        # update V
        Vn = (I + W) * V

        # calculate new average error
        for k in 1:K
            Cs[:, :, k] = Vn * Cs[:, :, k] * Vn'
            errs[s+1] += norm(Cs[:, :, k] - Diagonal(diag(Cs[:, :, k])))/K
        end

        if errs[s+1] > errs[s]
            break
        else
            V = Vn
        end
    end

    return V, errs[1:s+2]
end

# test fast_frobenius

# Cs = [rand(3, 3) for i in 1:10]
# V, errs = fast_frobenius(Cs, sweeps=100, theta=0.5)
# println(V)
# println(errs)

arr1 = [0 1 -2; 1 1 0; -2 0 3]
arr2 = [-1 3 -1; 3 5 -1; -1 -1 1]
arr = cat(arr1, arr2, dims=3)

V, errs = fast_frobenius(arr)
println(V)
# println(errs)
