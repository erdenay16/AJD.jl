module AJD


include("QDiag/_KN3_approach.jl")

import LinearAlgebra: isposdef, eigen, Diagonal, norm
import Random: Xoshiro

"""
    QDiag(C_0, C, weights, approach, tolerance, maximum_iteration, 
        random_number_generator = Xoshiro()) -> AbstractArray{<:Real}

This is an implementation of the algorithm introduced in: 
    Vollgraf, Roland & Obermayer, Klaus. (2006). 
    Quadratic optimization for simultaneous matrix diagonalization. 
    Signal Processing, IEEE Transactions on. 54. 3270 - 3278. 10.1109/TSP.2006.877673.
Namings of the parameters are consistent with the paper for easier understanding. 

# Arguments
- C_0::AbstractArray{<:Real}: This is the sphering matrix introduced as ``C^{(0)}``.
- C::AbstractArray{<:Real}: This is the set of matrices that is introduced as ``C``
- weights::AbstractArray{<:Real}:: This is the weights vector ``mathbf{\alpha}``. It will be
    later normalized such that ``sum{k=1}{K} \alpha_k = 1``.
- approach::String: This is the flag for approaches "NK3" and "N5" introduced in the paper.
- tolerance::Real: This is the tolerance for the error.
- maximum_iteration::Integer: Maximum number of iterations.
- random_number_generator::Xoshiro: This is random number generator for matrix ``W``. Default
    generator generates random numbers based on default_rng() but a seed can be introduced 
    by the user. This argument is added for testing purposes.
"""
function QDiag(
    C_0::AbstractArray{<:Real},
    C::AbstractArray{<:Real},
    weights::AbstractArray{<:Real},
    approach::String,
    tolerance::Real,
    maximum_iteration::Integer,
    random_number_generator::Xoshiro = Xoshiro(),
)::AbstractArray{<:Real}

    @assert isposdef(C_0) "C_0 must be positive-definite."
    @assert ((approach == "KN3") || (approach == "N5"))
    @assert length(size(C)) == 3 "C must be a tensor with dimensions K x M x N"
    @assert size(C_0, 1) == size(C_0, 2) == size(C, 1) == size(C, 2) "C must be a tensor with dimensions M x N x K and C_0 must be a matrix with dimensions M x N  where M = N"
    @assert length(size(weights)) == 1 "Weights variable must be a vector."
    @assert size(weights, 1) == size(C, 3) "Length of the weights must be equal to third dimension of C matrix."
    @assert tolerance >= 0 "Tolerance must be a nonnegative real number."
    @assert maximum_iteration > 0 "Maximum iteration must be positive."


    result = _optimize(
        C_0,
        C,
        weights,
        approach,
        tolerance,
        maximum_iteration,
        random_number_generator,
    )
end

function _optimize(
    C_0,
    C,
    weights,
    approach,
    tolerance,
    maximum_iteration,
    random_number_generator,
)
    if approach == "KN3"
        result = _apply_KN3(
            C_0,
            C,
            weights,
            tolerance,
            maximum_iteration,
            random_number_generator,
        )
    else
        # N5 will be implemented
    end
    return result
end

export QDiag

export ffdiag, plot_convergence, getW, off, cost_off, get_off, normit
end
