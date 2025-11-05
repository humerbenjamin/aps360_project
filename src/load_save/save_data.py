import os
import pickle

def save_single_equilibrium_state(data, statenumber):
    with open(f'data/rawdata/{statenumber}.pkl', 'wb') as f:
        pickle.dump(data, f)

    return

def initialize_data_saving_directories():
    if not os.path.isdir("data"):
        os.mkdir("data")
    if not os.path.isdir("data/rawdata"):
        os.mkdir("data/rawdata")

    print("Data Directory Structure Initialized")

    return