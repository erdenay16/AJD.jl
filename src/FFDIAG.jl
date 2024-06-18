
using LinearAlgebra: I, diag, diagm, norm

# Input: Set of matrices same size Cs
function ffdiag(Cs)
    
    # Init
    runs = 100
    # stepsize
    theta  = 0.5
    # turn set of mattrices into array
    Carr = collect(Cs)
    dim1, dim2 = size(Carr[1])
    
    # Cs must be set of square matrices -> check if square
    (dim1 == dim2) || throw(DimensionMismatch("Error: Matrices not square"))
    
    W = zeros(dim1, dim2)
    # V = I
    # init of V from joint diagonalizer python code
    V = [-0.85065081  0.52573111
     0.52573111  0.85065081]
    n = 1
    
    Ds = [zeros(dim1) for _ in Carr]
    Es = [zeros(dim1,dim2) for _ in Carr]

    errs = zeros(runs + 1)

    # doesnt do anything right now because V = I 
    for C in Carr
        C = V * C * V'
    end 


    # loop
    for run in 1:runs
        
        # fill the Ds and Es matrices
        for i in eachindex(Carr)
            
            Ds[i] = (diag(Carr[i]))
            Es[i] = copy(Carr[i])

            # make the diagonal of Es zero
            for j in 1:dim1
                Es[i][j,j] = 0
            end
        end
            
        # calculate W (17) 
        z = zeros(dim1, dim2)
        y = zeros(dim1, dim2)
        
        # Todo: Implement broadcast for the last loop sum(Ds[:,i] .* Ds[:,j] 
        for i in 1:dim1
            for j in 1:dim2
                for k in 1:length(Carr)
                    
                    D = Ds[k]
                    E = Es[k]

                    z[i, j] += (D[i] * D[j])
                    y[i, j] += 0.5 * D[j] * (E[i, j]+E[j, i]) # This is the same as += 0.5* Ds[k,j]*(Es[k,i,j] + Es[k,j,i]) as according to the paper

                end
            end
        end

        # Check if we need to the for loops are programmed correctly 
        # put the broadcast directly in the loop and elimnate the above 3 for loops
        for i in 1:dim1-1
            for j in (i+1):dim2
                if i != j
                    
                    W[i, j] = (z[i, j] * y[j, i] - z[i, i] * y[i, j]) / (z[j, j] * z[i, i] - z[i,j]^2)#(z[i, j]*y[i, j]-y[i, i]*z[i, j])/(z[j,j]*z[i,i]-z[i,j]^2)
                    W[j, i] = (z[i, j] * y[i, j] - z[j, j] * y[j, i]) / (z[j, j] * z[i, i] - z[i,j]^2)#(z[i,j]*y[i,j]-z[j,j]*y[j,i])/(z[j,j]*z[i,i]-z[i,j]^2)

                end
            end
        end


        # Compute the Frobenius norm of W
        # todo: Check why the python code uses a while loop here?
        if (norm(W)) > theta
            W = (theta/(norm(W)) * W)
        end

        V = (I + W) * V
    
        for i in eachindex(Carr)  
            Carr[i] = ((I + W) * Carr[i] * (I + W)') 
            
            errs[run + 1] += norm(Cs[i] - diagm(diag(Cs[i]))) / length(Carr)
        end  
        
        tol=1e-6
        # check for convergence
        if abs(errs[run + 1] - errs[run]) < tol
            @info "Converged at iteration ", run
            return Carr, V#errs[1:run + 1]
        end


        n += 1

 
        @info "C: ", Carr, V


    end


    @info "End: ", Carr, V


    return Carr, V

end


 D = [[1.0 2.0
 3.0 4.0], [4.0 5.0 
 5.0 6.0]]

ffdiag(D)