import os
import pickle

def save_single_equilibrium_state(data, statenumber, rawdata_folder="data/rawdata/"):
    with open(f'{rawdata_folder}{statenumber}.pkl', 'wb') as f:
        pickle.dump(data, f)

    return


def save_noisy_data_single_equilibrium_state(data, statenumber, noisedata_folder="data/noisedata/"):
    with open(f'{noisedata_folder}{statenumber}.pkl', 'wb') as f:
        pickle.dump(data, f)

    return


def initialize_data_saving_directories():
    if not os.path.isdir("data"):
        os.mkdir("data")
    if not os.path.isdir("data/rawdata"):
        os.mkdir("data/rawdata")
    if not os.path.isdir("data/noisedata"):
        os.mkdir("data/noisedata")

    print("Data Directory Structure Initialized")

    return