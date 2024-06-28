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
    # matrix = [matrixdepot("lehmer", 3) for i in 1:2]
    matrix = [matrixdepot("lehmer", 3), matrixdepot("lehmer", 3)]

    # plot(matrix[1])
    # matrix = matrixdepot("minij", 5)
    
    # matrix = [[0.0 1.0 -2.0
    #      1.0 1.0 0.0
    #      -2.0 0.0 3.0],
    #      [-1.0 3.0 -1.0
    #       3.0 5.0 -1.0
    #       -1.0 -1.0 1.0]]

    result = ffdiag(matrix)

    # println("The result is: $result")
    # @info(matrix[1])
    


    # is_diagonal = isdiag(result)

    # println("The result is a diagonal matrix: $is_diagonal")

end

# Test error, when matrices have different dimensions
# Test error, when matrices not square
# Test correct Ds is matrix with diagonal elements and and Es is matrix with all elements except diagonal is zero