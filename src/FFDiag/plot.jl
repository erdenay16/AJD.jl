using Plots


function plot_convergence(errs)

    plot(
    1:length(errs), errs;
    # Text labels
    title="Plotting Convergence",
    xlabel="number of iteration",
    ylabel="convergence error",
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