# import standard packages
    # import full packages
import sys
import os
import numpy as np
    # import specific functions
from scipy.stats import norm


# imports from within this project
    # properly defined functions
from src.load_save import load_single_equilibrium_state, load_noisedata_single_equilibrium_state, load_sensordata_single_equilibrium_state
from settings import load_synthetic_measurement_settings
from src.plotting import plot_F_p_over_psi_probability


def calculate_naive_bayes_probabilities_all_substates(rawdata_folder="data/rawdata/"):
    num_substates = load_synthetic_measurement_settings()["num_measurements_per_state"]
    
    # get list of rawdata files
    file_paths = []
    with os.scandir(rawdata_folder) as entries:
        for entry in entries:
            if entry.is_file():
                file_paths.append(entry.path)
    
    true_state_ranks = []
    for i in range(len(file_paths)):
        for j in range(num_substates):
            state_probabilities, state_numbers = calculate_naive_bayes_probabilities_single_substate(statenumber=i, substatenumber=j, file_paths=file_paths, supress_progress_bar=True)
            true_state_ranks.append(true_state_rank(state_probabilities, state_numbers, i))
        progress_bar(i+1, len(file_paths))
        i += 1

    print("\n\n")
    print(f"Number of Synthetic States Evaluated: {len(true_state_ranks)}")
    print(f"Number Correct: {true_state_ranks.count(0)}, Percent Correct: {(true_state_ranks.count(0)/len(true_state_ranks)):6.2f}%")
    print(f"Mean Distance From the Correct Value: {np.mean(true_state_ranks)}, Standard Deviation: {np.stdev(true_state_ranks)}")
    print(f"L2 Error: {np.mean(np.array(true_state_ranks) ** 2)}")



    return

def calculate_naive_bayes_probabilities_single_substate(statenumber=0, substatenumber=0, file_paths=[], rawdata_folder="data/rawdata/", printouts=False, plot_quantities=False, supress_progress_bar=False):
    noisy_data = load_noisedata_single_equilibrium_state(statenumber)[substatenumber]
    equil_state = load_single_equilibrium_state(0)
    if printouts:
        print("\nState Data:", equil_state['settings'], "\n")

    if len(file_paths) == 0:
        # get list of rawdata files
        file_paths = []
        with os.scandir(rawdata_folder) as entries:
            for entry in entries:
                if entry.is_file():
                    file_paths.append(entry.path)

    state_probabilities = []
    state_numbers = []
    for i in range(len(file_paths)):
        sensordata = load_sensordata_single_equilibrium_state(i)
        p = calculate_p_single_equilibrium(noisy_data, sensordata)

        state_probabilities.append(p)
        state_numbers.append(i)

        if not supress_progress_bar:
            progress_bar(i+1, len(file_paths), message="Naive Bayes")

    norm_coeff = np.nansum(state_probabilities)

    state_probabilities = state_probabilities / norm_coeff

    if printouts:
        print(f"\n{state_probabilities}")

    if plot_quantities:
        plot_F_p_over_psi_probability(state_probabilities, state_numbers, statenumber)

    return state_probabilities, state_numbers



def calculate_p_single_equilibrium(noisy_data, sensor_data):
    p=1
    for i in range(len(sensor_data)):
        for j in range(len(sensor_data[i])):
            offset = (noisy_data[i][j][1] - sensor_data[i][j]) / noisy_data[i][j][4]
            p = p*norm.pdf(offset, loc=0, scale=1)
            # print(p, offset, norm.pdf(offset, loc=0, scale=1))

    return p


def true_state_rank(state_probabilities, state_numbers, true_state):
    """
    Returns the rank (0-indexed) of the true_state after sorting 
    by descending state_probabilities.
    """
    # Sort indices by descending probability
    sorted_indices = np.argsort(state_probabilities)[::-1]
    
    # Sort the state_numbers accordingly
    sorted_states = np.array(state_numbers)[sorted_indices]
    
    # Find where true_state is in the sorted list
    rank = np.where(sorted_states == true_state)[0][0]
    
    return rank



# get the value of psi at a given R, Z pair for constants R0, a, b, c0
####################################################################################################
def progress_bar(progress, total, length=40, message=""):
    percent = 100 * (progress / total)
    filled = int(length * progress // total)
    bar = '█' * filled + '-' * (length - filled)
    sys.stdout.write(f'{message}\r|{bar}| {percent:6.2f}%')
    sys.stdout.flush()
    return
