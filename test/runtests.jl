using AJD
using Test

@testset "AJD.jl" begin
    include("QDiag_tests/QDiag_NK3_test.jl")
    include("FFDaig_tests/ffdiag_test_util.jl")
    include("FFDaig_tests/FFDiag_test.jl")
end
