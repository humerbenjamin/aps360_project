import matplotlib.pyplot as plt

def plot_psi_colormap(R_grid, Z_grid, psi_grid):
    plt.pcolormesh(R_grid, Z_grid, psi_grid, cmap="viridis_r")
    plt.colorbar(label=r"$\psi(R,Z)$")
    plt.xlabel("R")
    plt.ylabel("Z")
    plt.title(r"Poloidal Magnetic Flux $\psi(R,Z)$")
    plt.savefig("output_plots/psi_colormap.png")
    plt.close()

    return


def plot_psi_colormap_with_contours(R_grid, Z_grid, psi_grid, psi_contour_values, R_axis, statenumber=0):
    # plt.rcParams["font.family"] = "serif"
    plt.rcParams.update({
        "text.usetex": True,            # use LaTeX
        "font.family": "serif",
    })
    # set up figure
    plt.figure(figsize=(9,7), dpi=250)

    # plot quantities
    plt.vlines(x=R_grid[0][0], ymin=Z_grid[0][0], ymax=Z_grid[-1][0], color='k', linewidth=8, label="Inner Radial Boundary")
    plt.pcolormesh(R_grid, Z_grid, psi_grid, cmap="viridis_r")
    cbar = plt.colorbar()
    cbar.set_label(label=r"$\psi(R,Z)$ $[\mathrm{Wb/rad}]$", fontsize=24)
    cbar.ax.tick_params(labelsize=16)
    plt.contour(R_grid, Z_grid, psi_grid, levels=psi_contour_values, colors='r')
    plt.scatter(R_axis, 0, marker="x", color='r', label="Magnetic Axis")

    # set up axes
    plt.xlabel(r"R $[\mathrm{m}]$", fontsize=20)
    plt.ylabel(r"Z $[\mathrm{m}]$", fontsize=20)
    plt.xlim(xmin=R_grid[0][0], xmax=R_grid[-1][-1])
    plt.ylim(ymin=Z_grid[0][0], ymax=Z_grid[-1][-1])
    plt.tick_params(axis="both", which="major", labelsize=16)
    
    # title and legend
    plt.title(r"Poloidal Magnetic Flux $\psi(R,Z)$", fontsize=20)
    plt.plot([],[], color='r', label=r"Closed Contours of Constant $\psi$")
    plt.legend(fontsize=18)

    # layout and save figure
    plt.tight_layout()
    plt.savefig(f"output_plots/{statenumber}_psi_colormap_contours.png")
    plt.close()

    return