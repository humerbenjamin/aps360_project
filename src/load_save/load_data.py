import pickle

def load_single_equilibrium_state(statenumber):
    with open(f'data/rawdata/{statenumber}.pkl', 'rb') as f:
        data = pickle.load(f)

    return data