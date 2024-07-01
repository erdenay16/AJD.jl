using AJD: ffdiag, plot_convergence
using Test
using MatrixDepot: matrixdepot
#using Plots: plot
#using LinearAlgebra: I

include("ffdiag_test_util.jl")

@testset "AJD.jl" begin
    # Write your tests here.
    # generate 5 random matrices
    #=matrix = [matrixdepot("minij", 5) for i in 1:2]
    plot(matrix[1])
    # matrix = matrixdepot("minij", 5)

    result = ffdiag(matrix)

    @info(result[1])
    # @info(matrix[1])

    =#

    # is_diagonal = isdiag(result)

    # println("The result is a diagonal matrix: $is_diagonal")
    @test begin
        K = 15
        m_size = 10
        Cx, A = genFullyDiagMs(K, m_size)
        Cs, V, errs = ffdiag(Cx)
        check_diag(Cs)
    end

    @test begin
        K = 15
        m_size = 10
        Cx, Ax = genApproxDiagMs(K, m_size)
        Cs, V, errs = ffdiag(Cx)
        check_diag(Cs)
    end

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

        isapprox(V, res)
    end

    @test begin
        Ck = [1.0 2.0 3.0
            4.0 5.0 6.0
            7.0 8.0 9.0]
        Ck2 = [10.0 11.0 12.0
            13.0 14.0 15.0
            16.0 17.0 18.0]
        Cs = ones(3, 3, 2)
        Cs[:, :, 1] = Ck
        Cs[:, :, 2] = Ck2
        V = I
        cost_off(Cs, V)

        Vtest = [1.0 2.0 3.0
            4.0 5.0 6.0
            7.0 8.0 9.0]

        normit(Vtest)

        off(V, Ck)
        get_off(V, Cs)
    end
end

# Test error, when matrices have different dimensions
# Test error, when matrices not square
# Test correct Ds is matrix with diagonal elements and and Es is matrix with all elements except diagonal is zero
