using AJD: ffdiag, plot_convergence, getW, off, cost_off, get_off, off, normit
using LinearAlgebra: tr, I, opnorm
using Test
using MatrixDepot: matrixdepot
using Plots: plot


include("ffdiag_test_util.jl")


@testset "AJD.jl" begin
    
    @testset "test_helper_funktions" begin
        Cs = ones(2, 2, 2)
        Cs[:, :, 1] = [1.0 0.2; 0.2 0.8]
        Cs[:, :, 2] = [0.5 0.3; 0.3 0.5]
        V = [1.0 0.0; 0.0 1.0]
        # getW
        @test begin
            W = getW(Cs)
            res = [0.0 -2.0; 1.4 0.0]
            isapprox(W, res)
        end
        # off
        @test begin
            res = 0.2592
            isapprox(off(Cs[:, :, 1], V), res)
        end
        # get_off
        @test begin
            res = 0.2600
            isapprox(get_off(V, Cs), res)
        end
        # cost_off
        @test begin
            res = 0.2600
            isapprox(cost_off(Cs, V), res)
        end
        # normit
        @test begin
            V = [0.5 0.3; 0.3 0.5]
            res = [0.8575 0.5145
                0.5145 0.8575]
            normV = normit(V)
            isapprox(round.(normV, digits=4), res)
        end
    end

    @testset "ffdiag_assertions" begin
        Cs1 = ones(2, 2, 1)
        Cs1[:,:,1] = [1.0 0.2; 0.2 0.8]
        @test_throws AssertionError ffdiag(Cs1[:,:,1])
        @test_throws AssertionError ffdiag(Cs1)
        Cs2 = ones(2, 2, 2)
        Cs2[:, :, 1] = [1.0 0.2; 0.2 0.8]
        Cs2[:, :, 2] = [0.5 0.3; 0.3 0.5]
        @test_throws AssertionError ffdiag(Cs2, -10)
        @test_throws AssertionError ffdiag(Cs2, 100, -1.0)
        Cs3 = ones(2, 3, 2)
        Cs3[:, :, 1] = [1.0 0.2 0.2; 0.2 0.8 0.2]
        Cs3[:, :, 2] = [0.5 0.3 0.2; 0.3 0.5 0.2]
        @test_throws AssertionError ffdiag(Cs3)
        Cs4 = ones(2, 2, 2)
        Cs4[:, :, 1] = [1.0 0.2; 0.3 0.8]
        Cs4[:, :, 2] = [0.5 0.3; 0.3 0.5]
        @test_throws AssertionError ffdiag(Cs4)
    end


    @testset "test_diagonality" begin
        @test begin
            C0 = zeros(Float64, 5, 5, 2)

            C0[:, :, 1] = [
                2.235839 -0.132419 -0.468679 0.712025 -0.665429
                -0.132419 2.400627 -1.469896 -2.472759 0.231204
                -0.468679 -1.469896 3.090794 1.436710 0.265003
                0.712025 -2.472759 1.436710 2.876796 -0.093264
                -0.665429 0.231204 0.265003 -0.093264 1.118756
            ]

            C0[:, :, 2] = [
                -6.446048 -0.017503 3.762060 -1.269919 3.493108
                -0.017503 -0.481916 1.581665 0.790263 0.637676
                3.762060 1.581665 -3.591888 -0.743341 -1.778091
                -1.269919 0.790263 -0.743341 -1.202152 0.115708
                3.493108 0.637676 -1.778091 0.115708 -1.962549
            ]

            C, V, err = ffdiag(C0)

            res = [
                0.334704 0.538349 -0.177880 0.501605 -0.561164
                -0.204070 0.716255 0.037965 0.642288 -0.177085
                -0.194283 0.398450 0.442984 0.229722 -0.744638
                -0.423837 0.442813 -0.185234 0.675757 -0.365129
                0.465941 -0.626439 0.101427 -0.125307 0.603725
            ]

            isapprox(round.(V, digits=3), round.(res, digits=3))
        end

        # Makes 10 10x10 matrices and checks if the off-diagonal elements are zero
        @test begin
            K = 10
            n = 10

            C0, _ = genFullyDiagMs(n, K)
            _, V, _ = ffdiag(C0)

            f = 0

            for k in 1:K
                F = V * C0[:, :, k] * V'
                f += tr(F' * F) - tr(F .* F)
            end

            isapprox(0, f, atol=1e-9)
        end

        # Makes 10 10x10 matrices and checks if V* all C0 * V' are diagonal
        @testset "does_v_diganoalize_all" begin
            K = 10
            n = 10

            C0, _ = genFullyDiagMs(n, K)
            _, V, _ = ffdiag(C0)

            for k in 1:K
                F = V * C0[:, :, k] * V'
                @test isapprox(F, diagm(diag(F)), atol=1e-9)
            end
        end
    end

    @testset "plot_convergence_test" begin
        errs = [0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001]
        plot = plot_convergence(errs)
        @test !isnothing(plot)
    end
end



