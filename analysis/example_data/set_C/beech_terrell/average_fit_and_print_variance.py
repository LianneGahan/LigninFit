##### Added initial guess for test function parameters :partho
import numpy as np
import sys
import os
import json
from scipy import optimize
from latest.Code.compare_simulation_to_exp import calculate_cost
from latest.Code.figure_generation_functions import import_library, calculate_avg_bond_distribution

#argv[1] = file path to output
#argv[2] = Name of mean file
#argv[3] = name of file that contains the keywords (e.g. low,medium,high)

#simu_params = np.loadtxt("latest/Params/simulation_parameters.txt");
#N_files = int(simu_params[-3])
#print("N_files = " + str(N_files))
cost = 10

path = sys.argv[1] + "Output/library.json"
print(path)
if os.path.exists(path):
    # Import best fit data and calculate bond distribution
    library = import_library(path, "ligninkmc_")
    simulated_distribution = calculate_avg_bond_distribution(library[1])

    biomass_path = sys.argv[1] + "Params/biomass_data.json"
    # Open biomass JSON file for reading
    with open(biomass_path, 'r') as file:
        biomass_data = json.load(file)

        # Calculate cost between target and simulated bond distributions
        cost = calculate_cost(biomass_data, simulated_distribution)

print("var = " + str(cost))
with open(sys.argv[1] + "Output/var.txt","w") as f:
    f.write(str(cost))
    exit()
