# import standard packages
    # whole packages
import numpy as np
import sys
    # subpackages
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors

from src.load_save import load_single_equilibrium_state


def plot_F_p_over_psi(state_probabilities, state_numbers, true_state_number):
    plt.rcParams.update({
        "text.usetex": True,            # use LaTeX
        "font.family": "serif",
    })
    fig, axs = plt.subplots(1,2, figsize=(10,4), constrained_layout=True, dpi=300)

    # reorder state numbers
    indices = np.argsort(state_probabilities)
    state_probabilities = [state_probabilities[i] for i in indices]
    state_numbers = [state_numbers[i] for i in indices]

    # Set up colormap
    norm = mcolors.Normalize(vmin=np.min(state_probabilities), vmax=np.max(state_probabilities))
    cmap = cm.plasma  # choose any colormap
    colors = cmap(norm(state_probabilities))

    total_num_states = len(state_probabilities)

    for i in range(len(state_probabilities)):
        equilibrium = load_single_equilibrium_state(state_numbers[i])
        psi_1D = equilibrium['psi_1D']
        F_1D = equilibrium['F_1D']
        p_1D = equilibrium['p_1D']
        
        # F of psi
        axs[0].plot(psi_1D, F_1D, color=colors[i], linewidth=0.5)

        # P of psi
        axs[1].plot(psi_1D, p_1D, color=colors[i], linewidth=0.5)

        if state_numbers[i] == true_state_number:
            # F of psi
            axs[0].plot(psi_1D, F_1D, color='r', linestyle='--', label="True State", zorder=20)

            # P of psi
            axs[1].plot(psi_1D, p_1D, color='r', linestyle='--', label="True State", zorder=20)

        progress_bar(i+1, total_num_states, message="Plotting Probability Distribution of States in F and P")
            


    # set up plot labels and visuals
        # F
    axs[0].set_xlabel(r"$\psi\,[\mathrm{Wb}]$", fontsize=18)
    axs[0].set_ylabel(r"$F\,[\mathrm{Tm}]$", fontsize=18)
    axs[0].grid(True)
    axs[0].legend()
        # p
    axs[1].set_xlabel(r"$\psi\,[\mathrm{Wb}]$", fontsize=18)
    axs[1].set_ylabel(r"$p\,[\mathrm{Pa}]$", fontsize=18)
    axs[1].grid(True)
    axs[1].legend()
        # combined
    plt.suptitle(r"$F(\psi)$ and $p(\psi)$", fontsize=24)

    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, ax=axs, label=r'$P(C_k|\vec{x})$')

    plt.savefig(f"data/plots/state_prob_dist_F_p_psi.png")
    plt.close()


    return


# get the value of psi at a given R, Z pair for constants R0, a, b, c0
####################################################################################################
def progress_bar(progress, total, length=40, message=""):
    percent = 100 * (progress / total)
    filled = int(length * progress // total)
    bar = '█' * filled + '-' * (length - filled)
    sys.stdout.write(f'{message}\r|{bar}| {percent:6.2f}%')
    sys.stdout.flush()
    return