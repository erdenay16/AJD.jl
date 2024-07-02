using Plots: plot

"""
plot_convergence(errs) 

Plot the convergence error over iterations on a logarithmic scale.

# Arguments
- err::AbstractVector{<:AbstractFloat} An array with the errors, where the index corresponds to the number of iteration.

"""

function plot_convergence(errs::AbstractVector{<:Real})

    plot(
    1:length(errs), errs;
    # Text labels
    title="Plotting Convergence",
    xlabel="number of iteration",
    ylabel="convergence error",
    label="ffdiag",
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

