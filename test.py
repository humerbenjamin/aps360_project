import matplotlib.pyplot as plt

import src.solver.define_equilibria as define_equilibria
import src.solver.generate_synthetic_measurements as generate_synthetic_measurements
import src.models.naive_bayes.naive_bayes as naive_bayes



from src.load_save import initialize_data_saving_directories, clear_data_directories

import src.plotting.plot_sensors as plot_sensors




if __name__ == '__main__':
    # clear_data_directories()

    # define_equilibria.define_multiple_equilibria(plot_quantities=False)

    # generate_synthetic_measurements.generate_synthetic_measurements()

    # plot_sensors.plot_sensors_on_machine_geometry()

    # naive_bayes.calculate_naive_bayes_probabilities_single_substate(250, 0, printouts=True, plot_quantities=True)
    naive_bayes.calculate_naive_bayes_probabilities_all_substates()

