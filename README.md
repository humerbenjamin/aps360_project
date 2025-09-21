# aps360_project
Repository for the course project in APS360.

# ./src

All code (GS solver, and ML models).

## ./src/gs_solver

Grad Shafranov equilibrium solver.



# ./data

All generated data and models

## ./data/rawdata

Raw GS equilibrium outputs.

## ./data/input_data

GS data that has been processed into synthetic measurements.

### ./data/input_data/noiseless

GS data that has been processed into synthetic measurements directly from the GS equilibrium

### ./data/input_data/noisy

GS data that has been processed into synthetic measurements and had noise added to it (paired to the associated equilibrium so that outputs are appropriate).



# ./settings

Settings files for the various processes within the repo.