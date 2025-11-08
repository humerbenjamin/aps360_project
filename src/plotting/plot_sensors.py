import matplotlib.pyplot as plt


from settings.load_settings import load_machine_geometry, load_synthetic_measurement_settings


def plot_sensors_on_machine_geometry():
    machine_geometry = load_machine_geometry()
    synthetic_measurement_settings = load_synthetic_measurement_settings()

    plt.rcParams.update({
        "text.usetex": True,            # use LaTeX
        "font.family": "serif",
    })
    # set up figure
    plt.figure(figsize=(9,7), dpi=250)

    # plot quantities
    plt.scatter(synthetic_measurement_settings['T_sensors'][0], synthetic_measurement_settings['T_sensors'][1], s=300, marker="x", color='r', label=r"$T$ Measurement Locations", zorder=5)
    plt.scatter(synthetic_measurement_settings['n_sensors'][0], synthetic_measurement_settings['n_sensors'][1], s=300, marker='s', color='b', label=r'$n$ Measurement Locations', zorder=5)
    plt.scatter(synthetic_measurement_settings['b_sensors'][0], synthetic_measurement_settings['b_sensors'][1], s=300, marker='o', color='g', label=r'$B_{\theta}$ Measurement Locations', zorder=5)

    # set up axes
    plt.xlabel(r"R $[\mathrm{m}]$", fontsize=20)
    plt.ylabel(r"Z $[\mathrm{m}]$", fontsize=20)
    plt.xlim(xmin=machine_geometry['R_maj_min'], xmax=machine_geometry['R_maj_max'])
    plt.ylim(ymin=machine_geometry['Z_min'], ymax=machine_geometry['Z_max'])
    plt.tick_params(axis="both", which="major", labelsize=16)
    
    # title and legend
    plt.title(r"Sensor Locations in the Machine Geometry", fontsize=24)
    plt.legend(fontsize=18)

    # layout and save figure
    plt.tight_layout()
    plt.savefig(f"data/plots/sensors.png")
    plt.close()
    return