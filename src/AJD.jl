module AJD

include("QDiag/_KN3_approach.jl")
include("FFDiag/ffdiag.jl")

import Plots: plot
import LinearAlgebra: isposdef, eigen, Diagonal, norm, diag
import Random: Xoshiro

"""
    QDiag(C_0, C, weights, approach, tolerance, maximum_iteration, return_iteration_errors
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
- return_iteration_errors::Bool: If true error log wil be returned.
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
    return_iteration_errors::Bool = false,
    random_number_generator::Xoshiro = Xoshiro(),
)

    @assert isposdef(C_0) "C_0 must be positive-definite."
    @assert ((approach == "KN3") || (approach == "N5"))
    @assert length(size(C)) == 3 "C must be a tensor with dimensions K x M x N"
    @assert size(C_0, 1) == size(C_0, 2) == size(C, 1) == size(C, 2) "C must be a tensor with dimensions M x N x K and C_0 must be a matrix with dimensions M x N  where M = N"
    @assert length(size(weights)) == 1 "Weights variable must be a vector."
    @assert size(weights, 1) == size(C, 3) "Length of the weights must be equal to third dimension of C matrix."
    @assert tolerance >= 0 "Tolerance must be a nonnegative real number."
    @assert maximum_iteration > 0 "Maximum iteration must be positive."


    result, iteration_errors = _optimize(
        C_0,
        C,
        weights,
        approach,
        tolerance,
        maximum_iteration,
        random_number_generator,
    )
    if return_iteration_errors
        return result, iteration_errors
    else
        return result
    end
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
        result, iteration_errors = _apply_KN3(
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
    return result, iteration_errors
end

"""
plot_convergence(errs, label) 

Plot the convergence error over iterations on a logarithmic scale.

# Arguments
- errs::AbstractArray{<:Real} An array with the errors, where the index corresponds to the number of iteration.
- label::String Label of the figure.
"""
function plot_convergence(errs::AbstractArray{<:Real}, title::String)

    plot(
    1:length(errs), errs;
    # Text labels
    title=title,
    xlabel="number of iterations",
    ylabel="error",
    label = "error_log",
    # Line style
    color=:black,
    linewidth=2,
    # Axis settings
    yscale=:log10,
    # Other
    fontfamily="Helvetica",
    size=(600, 350),
    dpi=300,
)
end

export QDiag, ffdiag, plot_convergence, getW, off, cost_off, get_off, normit
end
