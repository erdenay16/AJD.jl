using AJD
using Test
using Random

import AJD: _find_P
import LinearAlgebra: Diagonal
import Distributions: Normal

const C0 = [1 0; 0 1]
const C = ones(2, 2, 2)
const weights = ones(2)
const approach = "KN3"
const tolerance = 0.0001
const maximum_iteration = 100
const return_iteration_errors = true

@testset "QDiag_NK3_test" begin
    @testset "QDiag_NK3_initialization" begin
        @testset "QDiag_NK3_positive_definite_assertion_test" begin
            C0_test = [1 2; 3 4]
            @test_throws AssertionError QDiag(
                C0_test,
                C,
                weights,
                approach,
                tolerance,
                maximum_iteration,
                return_iteration_errors
            )
        end

        @testset "QDiag_NK3_approach_assertion_test" begin
            approach_test = "ABC"
            @test_throws AssertionError QDiag(
                C0,
                C,
                weights,
                approach_test,
                tolerance,
                maximum_iteration,
                return_iteration_errors
            )
        end

        @testset "QDiag_NK3_size_length_assertion_test" begin
            C_test = ones(2, 2)
            @test_throws AssertionError QDiag(
                C0,
                C_test,
                weights,
                approach,
                tolerance,
                maximum_iteration,
                return_iteration_errors
            )
        end

        @testset "QDiag_NK3_square_matrix_assertion_test" begin
            C_test = ones(2, 1, 3)
            @test_throws AssertionError QDiag(
                C0,
                C_test,
                weights,
                approach,
                tolerance,
                maximum_iteration,
                return_iteration_errors
            )
        end

        @testset "QDiag_NK3_weigths_assertion_test" begin
            weights_test = ones(1)
            @test_throws AssertionError QDiag(
                C0,
                C,
                weights_test,
                approach,
                tolerance,
                maximum_iteration,
                return_iteration_errors
            )
            weights_test = ones(2, 2)
            @test_throws AssertionError QDiag(
                C0,
                C,
                weights_test,
                approach,
                tolerance,
                maximum_iteration,
                return_iteration_errors
            )
        end

        @testset "QDiag_NK3_tolerance_assertion_test" begin
            tolerance_test = -1
            @test_throws AssertionError QDiag(
                C0,
                C,
                weights,
                approach,
                tolerance_test,
                maximum_iteration,
                return_iteration_errors
            )
        end

        @testset "QDiag_NK3_maximum_iteration_assertion_test" begin
            maximum_iteration_test = -1
            @test_throws AssertionError QDiag(
                C0,
                C,
                weights,
                approach,
                tolerance,
                maximum_iteration_test,
                return_iteration_errors
            )
        end

        @testset "QDiag_NK3_find_p_test" begin
            identity = [1 0; 0 1]
            A = [2 1; 1 1.5]
            P_calculated = _find_P(A)
            @test isapprox(P_calculated * A * P_calculated', identity)
        end

        @testset "QDiag_NK3_approximately_diagonalizable_apply_test" begin
            Random.seed!(31)
            N = 4
            K = 10            
            weights_test = ones(K)
            random_vectors = randn(K, N)
            C0_test = Diagonal(ones(N))
            C_test = zeros(N,N,K)
            A  = randn(N,N)
            mean = 0.0
            stddev = 0.01

            for i in 1:K
                C_test[:, :, i] = A* Diagonal(random_vectors[i,:]) * A'
            end

            C_test_copy = deepcopy(C_test)
            desired_result =  [0.038843224659167785 0.40564274783350457 0.4203758223050812 -0.17178143487587333;
             0.5692856346947416 0.02761609491364972 -0.8514101835655733 -0.9329861572745857;
              0.09943954362553906 -0.5579816730150662 0.0754968054831004 0.11192040192145523; 
              0.815179027693811 0.7234277883874186 0.30444227647291466 -0.29580701918603725]

            random_number_generator = Xoshiro(7)
            result, iteration_errors = QDiag(
                C0_test,
                C_test,
                weights_test,
                approach,
                tolerance,
                maximum_iteration,
                return_iteration_errors,
                random_number_generator
            )

            @test isapprox(result, desired_result)
        end

        @testset "QDiag_NK3_approximately_diagonalizable_apply_test" begin
            Random.seed!(31)
            N = 4
            K = 10            
            weights_test = ones(K)
            random_vectors = randn(K, N)
            C0_test = Diagonal(ones(N))
            C_test = zeros(N,N,K)
            A  = randn(N,N)
            mean = 0.0
            stddev = 0.01

            for i in 1:K
                noise = rand(Normal(mean, stddev), size(A))
                C_test[:, :, i] = (A + noise) * Diagonal(random_vectors[i,:]) * (A + noise)'
            end

            C_test_copy = deepcopy(C_test)
            desired_result =  [0.04019543687211646 0.4060527565775023 0.4223882475955159 -0.17280885284805103; 
            0.5647693874738838 0.0323845918409765 -0.851185484846341 -0.9292730603642344; 
            0.1012979228592627 -0.5565385207653325 0.07579880761244097 0.10839105393435808; 
            0.8180211468241104 0.7241113671193407 0.30220188524332625 -0.30796113242685197]

            random_number_generator = Xoshiro(7)
            result, iteration_errors = QDiag(
                C0_test,
                C_test,
                weights_test,
                approach,
                tolerance,
                maximum_iteration,
                return_iteration_errors,
                random_number_generator
            )

            @test isapprox(result, desired_result)
        end
        @testset "QDiag_NK3_approximately_diagonalizable_apply_test" begin
            Random.seed!(31)
            N = 4
            K = 10            
            weights_test = ones(K)
            random_vectors = randn(K, N)
            C0_test = Diagonal(ones(N))
            C_test = zeros(N,N,K)
            A  = randn(N,N)
            mean = 0.0
            stddev = 0.01

            for i in 1:K
                noise = rand(Normal(mean, stddev), size(A))
                C_test[:, :, i] = (A + noise) * Diagonal(random_vectors[i,:]) * (A + noise)'
            end

            C_test_copy = deepcopy(C_test)
            desired_result =  [0.04019543687211646 0.4060527565775023 0.4223882475955159 -0.17280885284805103; 
            0.5647693874738838 0.0323845918409765 -0.851185484846341 -0.9292730603642344; 
            0.1012979228592627 -0.5565385207653325 0.07579880761244097 0.10839105393435808; 
            0.8180211468241104 0.7241113671193407 0.30220188524332625 -0.30796113242685197]

            random_number_generator = Xoshiro(7)
            result, iteration_errors = QDiag(
                C0_test,
                C_test,
                weights_test,
                approach,
                tolerance,
                maximum_iteration,
                return_iteration_errors,
                random_number_generator
            )

            @test isapprox(result, desired_result)
        end
    end
end
