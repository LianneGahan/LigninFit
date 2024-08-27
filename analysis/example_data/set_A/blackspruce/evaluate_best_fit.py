import numpy as np
import os
import json
from latest.Code.figure_generation_functions import import_data, calculate_cost, calculate_avg_bond_distribution



# Function to convert scientific notation to float
def sci_to_float(s):
    try:
        return float(s)
    except ValueError:
        # If the conversion fails, try to handle scientific notation
        parts = s.split("e")
        base = float(parts[0])
        exponent = float(parts[1].replace("+", ""))
        return base * 10**exponent


path = "BEST_FIT/best_Run/Output/library.json)"
if os.path.exists(path): # Read in best fit for library
    # Import best fit data and calculate bond distribution
    library = import_data(path, "ligninkmc_")
    simulated_distribution = calculate_avg_bond_distribution(library)

    # Open biomass JSON file for reading
    with open("BEST_FIT/best_Run/params/biomass_data.json", 'r') as file:
        biomass_data = json.load(file)

        # Calculate cost between target and simulated bond distributions
        cost = calculate_cost(biomass_data, simulated_distribution)


        # Print fit values
        print("-------------------------------")
        print(f"Euclidian Distance between target and simulated bond distributions: {cost:.4f}")