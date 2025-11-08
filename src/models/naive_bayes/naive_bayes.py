# import standard packages
    # import full packages
import os
import numpy as np
    # import specific functions
from scipy.stats import norm


# imports from within this project
    # properly defined functions
from src.load_save import load_single_equilibrium_state, load_noisedata_single_equilibrium_state, load_sensordata_single_equilibrium_state
from settings import load_synthetic_measurement_settings
    # misc functions







def calculate_naive_bayes_probabilities_single_substate(statenumber=0, substatenumber=0, rawdata_folder="data/rawdata/"):
    noisy_data = load_noisedata_single_equilibrium_state(statenumber)[substatenumber]

    # get list of rawdata files
    file_paths = []
    with os.scandir(rawdata_folder) as entries:
        for entry in entries:
            if entry.is_file():
                file_paths.append(entry.path)

    state_probabilities = []
    for i in range(len(file_paths)):
        sensordata = load_sensordata_single_equilibrium_state(i)
        p = calculate_p_single_equilibrium(noisy_data, sensordata)

        state_probabilities.append(p)

    print(state_probabilities)

    norm_coeff = np.nansum(state_probabilities)
    print(norm_coeff)

    state_probabilities = state_probabilities / norm_coeff

    print(state_probabilities)


    return state_probabilities




def calculate_p_single_equilibrium(noisy_data, sensor_data):
    p=1
    for i in range(len(sensor_data)):
        for j in range(len(sensor_data[i])):
            offset = (noisy_data[i][j][1] - sensor_data[i][j]) / noisy_data[i][j][4]
            p = p*norm.pdf(offset, loc=0, scale=1)
            # print(p, sensor_data[i][j], noisy_data[i][j])
            print(p, offset, norm.pdf(offset, loc=0, scale=1))

    print("\n\n")

    return p