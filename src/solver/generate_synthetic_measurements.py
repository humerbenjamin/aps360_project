# import standard packages
    # import full packages
import os
import numpy as np
    # import specific functions
from scipy.stats import norm


# imports from within this project
    # properly defined functions
from src.load_save import load_single_equilibrium_state, save_noisy_data_single_equilibrium_state
from settings import load_synthetic_measurement_settings
    # misc functions
from .define_equilibria import psi_of_R_Z


####################################################################################################
####################################################################################################
# MAIN FUNCTION TO GENERATE SYNTHETIC MEASUREMENTS
####################################################################################################
####################################################################################################

# main function to generate synthetic measurements
####################################################################################################
def generate_synthetic_measurements(rawdata_folder="data/rawdata/"):
    synthetic_measurement_settings = load_synthetic_measurement_settings()
    for i in range(len(synthetic_measurement_settings['n_sensors'][2])):
        synthetic_measurement_settings['n_sensors'][2][i] = float(synthetic_measurement_settings['n_sensors'][2][i])


    # get list of rawdata files
    file_paths = []
    with os.scandir(rawdata_folder) as entries:
        for entry in entries:
            if entry.is_file():
                file_paths.append(entry.path)
    
    for i in range(len(file_paths)):
        # load in equilibrium state
        equilibrium = load_single_equilibrium_state(i)

        # get values at different sensor locations
        T_values, n_values, b_values = get_values_at_sensor_positions(synthetic_measurement_settings, equilibrium)

        # generate synthetic measurments
        state_data_noisy = {}
        for n in range(synthetic_measurement_settings['num_measurements_per_state']):
            T_noisy, n_noisy, b_noisy = get_noisy_sensor_values(T_values, n_values, b_values, synthetic_measurement_settings)
            state_data_noisy[n] = [T_noisy, n_noisy, b_noisy]
        save_noisy_data_single_equilibrium_state(state_data_noisy, i)

    return 


####################################################################################################
####################################################################################################
# HELPER FUNCTIONS
####################################################################################################
####################################################################################################

# get values at the positions of sensors
####################################################################################################
def get_noisy_sensor_values(T_values, n_values, b_values, synthetic_measurement_settings):
    T_noisy = []
    n_noisy = []
    b_noisy = []

    mu, sigma = 0, 1


    # get temperature values
    for i in range(len(synthetic_measurement_settings['T_sensors'][0])):
        x = np.random.normal(mu, sigma)
        p = norm.pdf(x, mu, sigma)
        T_noisy.append([T_values[i] + x * synthetic_measurement_settings['T_sensors'][2][i], x, p, synthetic_measurement_settings['T_sensors'][2][i]])

    # get density values
    for i in range(len(synthetic_measurement_settings['n_sensors'][0])):
        x = np.random.normal(mu, sigma)
        p = norm.pdf(x, mu, sigma)
        n_noisy.append([n_values[i] + x * synthetic_measurement_settings['n_sensors'][2][i], x, p, synthetic_measurement_settings['n_sensors'][2][i]])


    # get values of the magnetic field
    for i in range(len(synthetic_measurement_settings['b_sensors'][0])):
        x = np.random.normal(mu, sigma)
        p = norm.pdf(x, mu, sigma)
        b_noisy.append([b_values[i] + x * synthetic_measurement_settings['b_sensors'][2][i], x, p, synthetic_measurement_settings['b_sensors'][2][i]])

    return T_noisy, n_noisy, b_noisy


# get values at the positions of sensors
####################################################################################################
def get_values_at_sensor_positions(synthetic_measurement_settings, equilibrium):
    # set up arrays for values
    T_values = []
    n_values = []
    b_values = []

    # get equilibrium settings
    equilibrium_settings = equilibrium['settings']

    # get temperature values
    for i in range(len(synthetic_measurement_settings['T_sensors'][0])):
        T_values.append(linearly_interpolate(psi_of_R_Z(synthetic_measurement_settings['T_sensors'][0][i], synthetic_measurement_settings['T_sensors'][1][i],
                            equilibrium_settings['R0'], equilibrium_settings['a'], equilibrium_settings['b'], equilibrium_settings['c0']),
                            equilibrium['psi_1D'], equilibrium['T_1D']))

    # get density values
    for i in range(len(synthetic_measurement_settings['n_sensors'][0])):
        n_values.append(linearly_interpolate(psi_of_R_Z(synthetic_measurement_settings['n_sensors'][0][i], synthetic_measurement_settings['n_sensors'][1][i],
                            equilibrium_settings['R0'], equilibrium_settings['a'], equilibrium_settings['b'], equilibrium_settings['c0']),
                            equilibrium['psi_1D'], equilibrium['n_1D']))

    # get values of the magnetic field
    for i in range(len(synthetic_measurement_settings['b_sensors'][0])):
        b_values.append(get_B_at_location(synthetic_measurement_settings['b_sensors'][0][i], synthetic_measurement_settings['b_sensors'][1][i],
                            equilibrium_settings['R0'], equilibrium_settings['a'], equilibrium_settings['b'], equilibrium_settings['c0']))

    return T_values, n_values, b_values


# get B at location
####################################################################################################
def get_B_at_location(R, Z, R0, a, b, c0):
    # calculate step size
    dR = 1e-3
    dZ = 1e-3

    # calculate derivatives
    dpsi_dR = (psi_of_R_Z(R+dR, Z, R0, a, b, c0) - psi_of_R_Z(R-dR, Z, R0, a, b, c0)) / (2*dR)
    dpsi_dZ = (psi_of_R_Z(R, Z+dZ, R0, a, b, c0) - psi_of_R_Z(R, Z-dZ, R0, a, b, c0)) / (2*dZ)

    # calculate B
    B = np.sqrt(dpsi_dR**2 + dpsi_dZ**2) / R

    return B


# linear interpolate
####################################################################################################
def linearly_interpolate(psi_value, psi_array, value_array):
    """
    Linearly interpolate to find the value corresponding to psi_value.
    Automatically converts input arrays to float arrays.
    """
    psi_array = np.asarray(psi_array, dtype=float)
    value_array = np.asarray(value_array, dtype=float)

    return float(np.interp(psi_value, psi_array, value_array))