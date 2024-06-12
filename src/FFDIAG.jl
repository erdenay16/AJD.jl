
using LinearAlgebra

# Input: Set of matrices same size Cs
function ffdiag(Cs)
    # C must be set of square matrices 
    # -> add check is square, check all matrices same dimensions
    
    # Init

    runs = 100
    theta  = 0.5
    Carr = collect(Cs)
    dim1, dim2 = size(Carr[1])
    if dim1 == dim2
        @info "dims matrices: ", dim1, dim2
    else 
        @error "matrix not square" 
    end
    
    W = zeros(dim1, dim2)
    V = I
    n = 1
    
    Ds = zeros(dim1, length(Carr))
    Es = zeros(dim1, dim2, length(Carr))

    # doiesnt do anything right now because V = I
    for C in Carr
        C = V * C * V'
    end 

    # loop
    for run in 1:runs
        
        # fill the Ds and Es matrices
        for k in eachindex(Carr)
            Ds[:,k] = diag(Carr[k])
            Es[:,:,k] = Carr[k]

            # make the diagonal of Es zero
            for j in 1:dim1
                Es[j,j,k] = 0                
                end
        
        # calculate W (17) 
        z = zeros(dim1, dim2)
        y = zeros(dim1, dim2)
        
        # Todo: Implement broadcast for the last loop sum(Ds[:,i] .* Ds[:,j] 
        for i in 1:dim1
            for j in 1:dim2
                for k in eachindex(Carr)
                    z[i, j] += Ds[i, k] * Ds[i, k]
                    y[i, j] += 0.5* Ds[j , k]*(Es[i,j,k] + Es[j,i,k]) 
                    # y[i, j] += Ds[j, k] * Es[i, j, k] # This is the same as += 0.5* Ds[k,j]*(Es[k,i,j] + Es[k,j,i]) as according to the paper
                end
            end
        end

        # Check if we need to the for loops are programmed correctyl 
        # put the broadcast directly in the loop and elimnate the above 3 for loops
        for i in 1:dim1
            for j in i:dim2
                # This was put in place so we dont divide by zero
                # Todo: Investigate why zero even shows up
                if z[j,j]*z[i,i]-z[i,j]^2 != 0
                    W[i, j] = (z[i, j]*y[j, i]-y[i, i]*z[i, j])/(z[j,j]*z[i,i]-z[i,j]^2)
                    W[j, i] = (z[i,j]*y[i,j]-z[j,j]*y[j,i])/(z[j,j]*z[i,i]-z[i,j]^2)
                else
                    W[i, j] = 0
                    W[j, i] = 0
                end
            end
        end

        # Compute the Frobenius norm of W
        # todo: Check why the python code uses a while loop here?

        if (norm(W)) > theta
            W = (theta/(norm(W)) * W)
        end

        V = (I + W) * V

        @info(W)
        for K in eachindex(Carr)
            Carr[K] = (I + W[k]) * Carr[k] * (I + W[k])'
        end 

        n += 1

    end

    return Carr, V
    
end
end

export ffdiag