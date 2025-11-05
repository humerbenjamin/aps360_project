import os

from src.load_save import load_single_equilibrium_state
from settings import load_synthetic_measurement_settings



def generate_synthetic_measurements(rawdata_folder="data/rawdata/"):

    # get list of rawdata files
    file_paths = []
    with os.scandir(rawdata_folder) as entries:
        for entry in entries:
            if entry.is_file():
                file_paths.append(entry.path)
    
    print(file_paths)

    return



