import matplotlib.pyplot as plt

import src.solver.define_equilibria as define_equilibria
from src.load_save import initialize_data_saving_directories



if __name__ == '__main__':
    # initialize_data_saving_directories()

    define_equilibria.define_multiple_equilibria(True)

