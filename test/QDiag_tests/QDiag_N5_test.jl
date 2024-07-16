using AJD
using Test
using LinearAlgebra
using Random


function generate_pos_def_matrix(N, rng)
    A = randn(rng, N, N)
    return A' * A + I * N
end


@testset "Basic Functionality" begin
    N = 3
    K = 2
    rng = MersenneTwister(123)
    C_0 = generate_pos_def_matrix(N, rng)
    C = [generate_pos_def_matrix(N, rng) for _ in 1:K]
    C = reshape(C, N, N, K)
    weights = rand(rng, K)
    tolerance = 1e-6
    max_iter = 1000

    W = QDiag_N5(C_0, C, weights, tolerance, max_iter, rng)
    
    @test size(W) == (N, N)
end


@testset "Convergence with Tolerance" begin
    N = 3
    K = 3
    rng = MersenneTwister(123)
    C_0 = generate_pos_def_matrix(N, rng)
    C = [generate_pos_def_matrix(N, rng) for _ in 1:K]
    C = reshape(C, N, N, K)
    weights = rand(rng, K)
    tolerance = 1e-8
    max_iter = 1000

    W = QDiag_N5(C_0, C, weights, tolerance, max_iter, rng)
    
    @test size(W) == (N, N)
end


@testset "Maximum Iterations" begin
    N = 4
    K = 2
    rng = MersenneTwister(123)
    C_0 = generate_pos_def_matrix(N, rng)
    C = [generate_pos_def_matrix(N, rng) for _ in 1:K]
    C = reshape(C, N, N, K)
    weights = rand(rng, K)
    tolerance = 1e-12
    max_iter = 5

    W = QDiag_N5(C_0, C, weights, tolerance, max_iter, rng)
    
    @test size(W) == (N, N)
end


@testset "Weights Normalization" begin
    N = 3
    K = 3
    rng = MersenneTwister(123)
    C_0 = generate_pos_def_matrix(N, rng)
    C = [generate_pos_def_matrix(N, rng) for _ in 1:K]
    C = reshape(C, N, N, K)
    weights = [1.0, 2.0, 3.0]
    tolerance = 1e-6
    max_iter = 1000

    W = QDiag_N5(C_0, C, weights, tolerance, max_iter, rng)
    
    @test size(W) == (N, N)
    @test abs(sum(weights) - 1.0) < 1e-12
end
