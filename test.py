import matplotlib.pyplot as plt

import src.solver.define_equilibria as define_equilibria
import src.solver.generate_synthetic_measurements as generate_synthetic_measurements
from src.load_save import initialize_data_saving_directories




if __name__ == '__main__':
    # initialize_data_saving_directories()

    # define_equilibria.define_multiple_equilibria(True)

    generate_synthetic_measurements.generate_synthetic_measurements()