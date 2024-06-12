using AJD
using Test
using MatrixDepot
# include("../src/FFDIAG.jl")
# using .FFDIAG: ffdiag

import AJD: ffdiag

@testset "AJD.jl" begin
    # Write your tests here.
    # generate 5 random matrices
    matrix = [matrixdepot("minij", 5) for _ in 1:2]
    # matrix = matrixdepot("minij", 5)s<
    
    result = ffdiag(matrix)
    is_diagonal = isdiag(result)

    println("The result is a diagonal matrix: $is_diagonal")

end
