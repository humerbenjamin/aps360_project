import matplotlib.pyplot as plt


def plot_P_F_over_psi(psi_1D, F_1D, p_1D, statenumber=0):
    fig, axs = plt.subplots(1,2, figsize=(10,4))
    
    # F of psi
    axs[0].plot(psi_1D, F_1D)
    axs[0].set_xlabel(r"$\psi\,[\mathrm{Wb}]$", fontsize=18)
    axs[0].set_ylabel(r"$F\,[\mathrm{Tm}]$", fontsize=18)
    axs[0].grid(True)

    # P of psi
    axs[1].plot(psi_1D, p_1D)
    axs[1].set_xlabel(r"$\psi\,[\mathrm{Wb}]$", fontsize=18)
    axs[1].set_ylabel(r"$p\,[\mathrm{Pa}]$", fontsize=18)
    axs[1].grid(True)

    plt.suptitle(r"$F(\psi)$ and $p(\psi)$", fontsize=24)
    plt.tight_layout()
    plt.savefig(f"output_plots/{statenumber}_F_p_psi.png")
    plt.close()


    return