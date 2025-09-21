import yaml

def load_machine_geometry(printout=False):
    with open("settings/machine_geometry.yml", "r") as f:
        data = yaml.safe_load(f)
    if printout:
        print(data)

    return data

def load_solver_parameters(printout=False):
    with open("settings/solver_parameters.yml", "r") as f:
        data = yaml.safe_load(f)
    if printout:
        print(data)
    
    return data

def load_synthetic_measurements(printout=False):
    with open("settings/synthetic_measurements.yml", "r") as f:
        data = yaml.safe_load(f)
    if printout:
        print(data)
    
    return data


if __name__ == '__main__':
    load_machine_geometry(printout=True)
    load_solver_parameters(printout=True)
    load_synthetic_measurements(printout=True)