import os
import numpy as np

from src.load_save import load_single_equilibrium_state
from settings import load_synthetic_measurement_settings

from .define_equilibria import psi_of_R_Z



def generate_synthetic_measurements(rawdata_folder="data/rawdata/"):
    synthetic_measurement_settings = load_synthetic_measurement_settings()
    print(synthetic_measurement_settings)

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

        

        
    print(equilibrium)


    print(T_values, n_values, b_values)

    return 


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


def linearly_interpolate(psi_value, psi_array, value_array):
    """
    Linearly interpolate to find the value corresponding to psi_value.
    Automatically converts input arrays to float arrays.
    """
    psi_array = np.asarray(psi_array, dtype=float)
    value_array = np.asarray(value_array, dtype=float)

    return float(np.interp(psi_value, psi_array, value_array))