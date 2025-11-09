import matplotlib.pyplot as plt


def plot_P_F_over_psi(psi_1D, F_1D, p_1D, statenumber=0):
    plt.rcParams.update({
        "text.usetex": True,            # use LaTeX
        "font.family": "serif",
    })
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
    plt.savefig(f"data/plots/{statenumber}_F_p_psi.png")
    plt.close()

    return


def plot_n_T_over_psi(psi_1D, n_1D, T_1D, statenumber=0):
    plt.rcParams.update({
        "text.usetex": True,            # use LaTeX
        "font.family": "serif",
    })
    fig, axs = plt.subplots(1,2, figsize=(10,4))
    
    # F of psi
    axs[0].plot(psi_1D, n_1D)
    axs[0].set_xlabel(r"$\psi\,[\mathrm{Wb}]$", fontsize=18)
    axs[0].set_ylabel(r"$n\,[\mathrm{m^{-3}}]$", fontsize=18)
    axs[0].grid(True)

    # P of psi
    axs[1].plot(psi_1D, T_1D)
    axs[1].set_xlabel(r"$\psi\,[\mathrm{Wb}]$", fontsize=18)
    axs[1].set_ylabel(r"$T\,[\mathrm{eV}]$", fontsize=18)
    axs[1].grid(True)

    plt.suptitle(r"$n(\psi)$ and $T(\psi)$", fontsize=24)
    plt.tight_layout()
    plt.savefig(f"data/plots/{statenumber}_n_T_psi.png")
    plt.close()

    return