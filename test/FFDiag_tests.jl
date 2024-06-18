using AJD
using Test
using MatrixDepot: matrixdepot
using Plots: plot
using LinearAlgebra: I
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

# Test error, when matrices have different dimensions
# Test error, when matrices not square
# Test correct Ds is matrix with diagonal elements and and Es is matrix with all elements except diagonal is zero