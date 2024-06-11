
using LinearAlgebra

# Input: Set of matrices same size Cs
function ffdiag(Cs)
    # C must be set of square matrices 
    # -> add check is square, check all matrices same dimensions
    
    # Init

    runs = 100
    theta  = 0.5
    Carr = collect(Cs)
    dim1,dim2 = size(Carr[1])
    if dim1 == dim2
        @info "dims matrices: ", dim1, dim2
    else 
        @error "matrix not square" 
    end
    
    W = zeros(dim1, dim2)
    V = I
    n = 1
    
    for C in Carr
        C = V * C * transpose(V)
    end 

    # loop
    for run in 1:runs
        
        # calculate W (17) or (18)


        if (norm(W)) > theta
            W = (theta/(norm(W)) * W)
        end

        V = (I + W) * V
        for C in Carr
            C = (I + W) * C * transpose((I + W))
        end 
        n += 1

    end

    return Carr, V

end

ffdiag(Cs)[]
