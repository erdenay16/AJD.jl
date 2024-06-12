using AJD
using Test
using MatrixDepot
using Plots
# include("../src/FFDIAG.jl")
# using .FFDIAG: ffdiag

import AJD: ffdiag

@testset "AJD.jl" begin
    # Write your tests here.
    # generate 5 random matrices
    matrix = [matrixdepot("minij", 5) for i in 1:2]
    plot(matrix[1])
    # matrix = matrixdepot("minij", 5)
    
    result = ffdiag(matrix)

    @info(result[1])
    # @info(matrix[1])
    


    # is_diagonal = isdiag(result)

    # println("The result is a diagonal matrix: $is_diagonal")

end
