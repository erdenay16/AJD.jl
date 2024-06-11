module AJD
using LinearAlgebra
# Write your package code here.

function FFDiag(cᵏ)

    sweeps = 100
    theta = 0.5

    # INPUT: Ck { Matrices to be diagonalized}
    K, m, N = size(cᵏ)
    W = zeros(K, m, m)

    # We are clever guessing the initial value of V
    _, V = eigen(Cs[1])

    z = zeros(K, m, m)
    y = zeros(K, m, m)

    
    Ds = zeros(m, m)
    Es = zeros(k, m, m)
    
    errors = zeros(K) 
    #initali

    for s in 1:sweeps 

        z = zeros(m, m)
        y = zeros(m, m)
        
        # Todo: Implement broadcast for the last loop sum(Ds[:,i] .* Ds[:,j] 
        for i in 1:m
            for j in 1:m
                for k in 1:K
                    z[i, j] += Ds[k, i] * Ds[k, j]
                    y[i, j] += Ds[k, j] * Es[k, i, j] # This is the same as += 0.5* Ds[k,j]*(Es[k,i,j] + Es[k,j,i]) as according to the paper
                end
            end
        end


        # Check if we need to the for loops are programmed correctyl 
        # put the broadcast directly in the loop and elimnate the above 3 for loops
        for i in 1:m-1
            for j in i+1:m
                W[i, j] = (z[i, j]*y[i, j]-y[i, i]*z[i, j])/(z[j,j]*z[i,i]-z[i,j]^2)
                W[j, i] = (z[i,j]*y[i,j]-z[j,j]*y[j,i])/(z[j,j]*z[i,i]-z[i,j]^2)
            end
        end

        # Compute the Frobenius norm of W
        # todo: Check why the python code uses a while loop here?
        Wf = norm(W, "fro")
        if Wf > theta
            W = theta * W / Wf
        end

        # Update V
        V = (I + W) * V

        for k in 1:K
            cᵏ[k] = (I + W) * cᵏ[k] * (I + W)'
            errors[s+1] += norm(cᵏ[k] - V * Ds * V', "fro")
        end

        if errors[s+1] > errors[s]
            break
        end
    end
end

export FFDiag
end
