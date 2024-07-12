using AJD
using Test

import AJD: _find_P
import Random: Xoshiro

const C0 = [1 0; 0 1]
const C = ones(2, 2, 2)
const weights = ones(2)
const approach = "KN3"
const tolerance = 0.0001
const maximum_iteration = 100

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
            )
            weights_test = ones(2, 2)
            @test_throws AssertionError QDiag(
                C0,
                C,
                weights_test,
                approach,
                tolerance,
                maximum_iteration,
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
            )
        end

        @testset "QDiag_NK3_find_p_test" begin
            identity = [1 0; 0 1]
            A = [2 1; 1 1.5]
            P_calculated = _find_P(A)
            @test isapprox(P_calculated' * A * P_calculated, identity)
        end

        @testset "QDiag_NK3_apply_test" begin
            C0_test = [20 12; 12 8]
            C_test = ones(2, 2, 3)
            weights_test = ones(3)
            desired_result = [
                1.0 0.8607279638770475
                0.8607279638770475 1.0
            ]
            C_test[:, :, 1] = [5 4; 4 5]
            C_test[:, :, 2] = [13 12; 12 13]
            C_test[:, :, 3] = [5 6; 6 8]

            random_number_generator = Xoshiro(7)
            return_iteration_errors = true
            result, iteration_errors = QDiag(
                C0_test,
                C_test,
                weights_test,
                approach,
                tolerance,
                maximum_iteration,
                random_number_generator,
                return_iteration_errors
            )

            @test isapprox(result * C0_test * result', desired_result)
            @test isequal(length(iteration_errors), 46)
        end
    end
end
