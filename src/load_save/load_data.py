import pickle

def load_single_equilibrium_state(statenumber, state_folder_path='data/rawdata/'):
    with open(f'{state_folder_path}{statenumber}.pkl', 'rb') as f:
        data = pickle.load(f)

    return data

def load_noisedata_single_equilibrium_state(statenumber, state_folder_path='data/noisedata/'):
    with open(f'{state_folder_path}{statenumber}.pkl', 'rb') as f:
        data = pickle.load(f)

    return data