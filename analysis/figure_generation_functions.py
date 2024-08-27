import numpy as np
import json
from scipy.optimize import curve_fit
#from rdkit import AllChem, Descriptors#
#import rdkit

#from rdkit import Chem
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import rdMolDescriptors
#from rdkit.Chem.AllChem import GetRDKitFPGenerator

from scipy import stats

def load_biomass_data(file_path):
    with open(file_path, 'r') as f:
        data = json.load(f)
    return data


def calculate_single_fingerprint(data):
    molecule = AllChem.MolFromSmiles(data["smilestring"])
    fingerprint = AllChem.GetMorganFingerprint(molecule, radius=2)
    return fingerprint

def calculate_fingerprints(lignin_library):
    sorted_data = sorted(lignin_library, key=lambda x: x['DP'])

    lignin_library = sorted_data[:2000]

    #print(lignin_library)
   # Generate fingerprint using an alternative method
    #print(rdkit.__version__)
    i=0

    fingerprints = []
    for data in lignin_library:
        if data["DP"]>2:
            print(i)
            molecule = AllChem.MolFromSmiles(data["smilestring"])
            fingerprints.append(AllChem.GetMorganFingerprint(molecule, radius=2))
    #molecules = [AllChem.MolFromSmiles(data["smilestring"]) for data in lignin_library]
    #fingerprints = [AllChem.GetMorganFingerprint(x) for x in molecules]
            i+=1
    return(fingerprints)

def compute_similarity_matrix(fps):
    n = len(fps)
    sim_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
            sim_matrix[i, j] = sim
            sim_matrix[j, i] = sim
    return sim_matrix

def compute_similarity_matrix_twolibraries(fps1, fps2):
    n1 = len(fps1)
    n2 = len(fps2)
    sim_matrix = np.zeros((n1, n2))
    for i in range(n1):
        for j in range(n2):
            sim = DataStructs.TanimotoSimilarity(fps1[i], fps2[j])
            sim_matrix[i, j] = sim
            #sim_matrix[j, i] = sim
    return sim_matrix

def compute_average_similarity(fps):
    n = len(fps)
    sims = []
    for i in range(n):
        for j in range(i, n):
            sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
            sims.append(sim)
    avg = np.sum(np.array(sims))/len(sims)#
    return(avg)

def euclidian_dist(biomass, simulated_distribution):
    """ Calculates the Euclidian distance between simulated and desired distribution 
    
    Input: 

    Biomass - a JSON file containing a dictionary with information about the biomass. This will
                contain a bond distribution, sg ratio, average molecular weight and the level of
                detail in the biomass data. i.e. Complete, Klose or Terrell type

    simulated distribution - A python dictionary containing propotions of bo4, bb, b5, b1, ao4, 4o5 and 55
                bond types


    Output:

    Cost - The euclidian distance between the simulated and target bond distribution
    """

    if biomass["bond_detail_level"] == "klose":
        _target_distribution = np.array([biomass["bo4"], 
                                        biomass["bb"],
                                        biomass["b5"]])
        _simulated_distribution = np.array([simulated_distribution["bo4"], 
                                    simulated_distribution["bb"],
                                    simulated_distribution["b5"]])
        
        # Normalise so we keep the ratios of bond types
        _target_distribution = _target_distribution/np.sum(_target_distribution)
        _simulated_distribution = _simulated_distribution/np.sum(_simulated_distribution)
        
    elif biomass["bond_detail_level"] == "Complete":
        # Calculate the difference between simulated and desired distribution
        _target_distribution = np.array([biomass["bo4"], 
                                        biomass["bb"],
                                        biomass["b5"],
                                        biomass["b1"],
                                        biomass["4o5"], 
                                        biomass["ao4"],
                                        biomass["55"]])
        
        _simulated_distribution = np.array([simulated_distribution["bo4"], 
                                    simulated_distribution["bb"],
                                    simulated_distribution["b5"],
                                    simulated_distribution["b1"],
                                   simulated_distribution["5o4"],
                                    simulated_distribution["ao4"], 
                                    simulated_distribution["55"]])
        
        # Normalise so we keep the ratios of bond types
        _target_distribution = _target_distribution/np.sum(_target_distribution)
        _simulated_distribution = _simulated_distribution/np.sum(_simulated_distribution)

    elif biomass["bond_detail_level"] == "terrell":
        # Calculate the difference between simulated and desired distribution
        _target_distribution = np.array([biomass["bo4"], 
                                        biomass["bb"],
                                        biomass["b5"],
                                        biomass["b1"],
                                        biomass["4o5"] + 
                                        biomass["ao4"] +
                                        biomass["55"]])
        
        _simulated_distribution = np.array([simulated_distribution["bo4"], 
                                    simulated_distribution["bb"],
                                    simulated_distribution["b5"],
                                    simulated_distribution["b1"],
                                    simulated_distribution["5o4"] +
                                    simulated_distribution["ao4"] +
                                    simulated_distribution["55"]])
        
        # Normalise so we keep the ratios of bond types
        _target_distribution = _target_distribution/np.sum(_target_distribution)
        _simulated_distribution = _simulated_distribution/np.sum(_simulated_distribution)

    euclidian_distance = np.linalg.norm(_target_distribution - _simulated_distribution)

    return(euclidian_distance)

def bond_metrics(biomass, simulated_distribution):
    """ Calculates the Euclidian distance between simulated and desired distribution 
    
    Input: 

    Biomass - a JSON file containing a dictionary with information about the biomass. This will
                contain a bond distribution, sg ratio, average molecular weight and the level of
                detail in the biomass data. i.e. Complete, Klose or Terrell type

    simulated distribution - A python dictionary containing propotions of bo4, bb, b5, b1, ao4, 4o5 and 55
                bond types


    Output:

    Cost - The euclidian distance between the simulated and target bond distribution
    """

    if biomass["bond_detail_level"] == "klose":
        _target_distribution = np.array([biomass["bo4"], 
                                        biomass["bb"],
                                        biomass["b5"]])
        _simulated_distribution = np.array([simulated_distribution["bo4"], 
                                    simulated_distribution["bb"],
                                    simulated_distribution["b5"]])
        
        # Normalise so we keep the ratios of bond types
        _target_distribution = _target_distribution/np.sum(_target_distribution)
        _simulated_distribution = _simulated_distribution/np.sum(_simulated_distribution)
        
    elif biomass["bond_detail_level"] == "Complete":
        # Calculate the difference between simulated and desired distribution
        _target_distribution = np.array([biomass["bo4"], 
                                        biomass["bb"],
                                        biomass["b5"],
                                        biomass["b1"],
                                        biomass["4o5"], 
                                        biomass["ao4"],
                                        biomass["55"]])
        
        _simulated_distribution = np.array([simulated_distribution["bo4"], 
                                    simulated_distribution["bb"],
                                    simulated_distribution["b5"],
                                    simulated_distribution["b1"],
                                   simulated_distribution["5o4"],
                                    simulated_distribution["ao4"], 
                                    simulated_distribution["55"]])
        
        # Normalise so we keep the ratios of bond types
        _target_distribution = _target_distribution/np.sum(_target_distribution)
        _simulated_distribution = _simulated_distribution/np.sum(_simulated_distribution)

    elif biomass["bond_detail_level"] == "terrell":
        # Calculate the difference between simulated and desired distribution
        _target_distribution = np.array([biomass["bo4"], 
                                        biomass["bb"],
                                        biomass["b5"],
                                        biomass["b1"],
                                        biomass["4o5"] + 
                                        biomass["ao4"] +
                                        biomass["55"]])
        
        _simulated_distribution = np.array([simulated_distribution["bo4"], 
                                    simulated_distribution["bb"],
                                    simulated_distribution["b5"],
                                    simulated_distribution["b1"],
                                    simulated_distribution["5o4"] +
                                    simulated_distribution["ao4"] +
                                    simulated_distribution["55"]])
        
        # Normalise so we keep the ratios of bond types
        _target_distribution = _target_distribution/np.sum(_target_distribution)
        _simulated_distribution = _simulated_distribution/np.sum(_simulated_distribution)

    euclidian_distance = np.linalg.norm(_target_distribution - _simulated_distribution)
    spearman_correlation = stats.spearmanr(_target_distribution, _simulated_distribution)
    pearson_correlation = stats.pearsonr(_target_distribution, _simulated_distribution)


    coefficient_dictionary = {
        "euclidian": euclidian_distance,
        "spearman": spearman_correlation.statistic,
        "pearson": pearson_correlation.statistic
    }
    return(coefficient_dictionary)

def import_library(json_file_name, string):
    """
    Input: File location of .json file containing the outputs of a Lignin KMC simulation
    Output: 
        _ids: The list of id keys associated to each dictionary (one lignin molecule)
        _dict: The library containing multiple dictionaries, each of which corresponds to the properties of one lignin molecule.
        Each dict contains:

        SMILES string of the molecule
        Bonds :"bo4" "bb" "b5" "b1" "a04" "55" "4o5" contributions
        func_grps
        MW
        DP

    """
# Create a defaultdict with a list as the default value
    _dict = []
    _ids = []
    # Load the JSON data
    with open(json_file_name, 'r') as json_file:
        json_data = json.load(json_file)

        for i in range(1000000):
            id_string = string + str(i)
            if id_string in json_data.keys():
                _dict.append(json_data[id_string])
                _ids.append(id_string)
            else: 
                if i!=0:
                    break
    return(_ids, _dict)

def calculate_avg_bond_distribution(library):
    """
    Calculate the normalised distribution of specific bond types within a collection of simulated structures.

    Inputs:
        simulation_result (list): A collection of simulation data containing adjacency matrices and bond information.
        num_sims (int): The number of simulations or structures in the provided data.

    Returns:
        dict: A dictionary containing the normalized distribution of bond types within the simulated structures.
            The keys are bond type names, and the values are the normalized counts of each bond type."""
    #get all bond dictionaries
    Bonds = [inner_dict["Bonds"] for inner_dict in library if "Bonds" in inner_dict]
    #print(Bonds)
    BO4 = []
    BB = []
    B5 = []
    _55 = []
    _4O5 = []
    B1 = []
    AO4 =[]
    
    # get all bond types
    for dictionary in Bonds:
        if 'bo4' in dictionary:
            BO4.append(dictionary['bo4'])
        if 'bb' in dictionary:
            BB.append(dictionary['bb'])
        if 'b5' in dictionary:
            B5.append(dictionary['b5'])
        if '55' in dictionary:
            _55.append(dictionary['55'])
        if '5o4' in dictionary:
            _4O5.append(dictionary['5o4'])
        if 'b1' in dictionary:
            B1.append(dictionary['b1'])
        if 'ao4' in dictionary:
            AO4.append(dictionary['ao4'])
    #Normalise the distribution:
    norm = np.sum(BO4) + np.sum(BB) + np.sum(B5) + np.sum(_55) + np.sum(AO4) + np.sum(_4O5) + np.sum(B1)

    dict = {
            "bo4": np.sum(BO4)/norm,
            "bb": np.sum(BB)/norm,
            "b5": np.sum(B5)/norm,
            "ao4": np.sum(AO4)/norm,
            "5o4": np.sum(_4O5)/norm,
            "b1": np.sum(B1)/norm, 
            "55": np.sum(_55)/norm
            }
    return(dict)

def calculate_bond_distribution_fitness(simulated_distribution, target_distribution):
    """
    Inputs:
    simulated_distribution - a dictionary which contains a value for each of the dominant bonds from a generated library of lignin. 
                            This does not have to be normalised but it is better if it is. 

    target_distribution - a dictionary which contains a value for each of the dominant bonds. This does not have to be 
                            normalised but it is better if it is.
    """
    _target_distribution = np.array([target_distribution["bo4"], 
                                    target_distribution["bb"],
                                    target_distribution["b5"],
                                    target_distribution["b1"],
                                    target_distribution["4o5"],
                                    target_distribution["ao4"],
                                    target_distribution["55"]])
    _target_distribution = _target_distribution/np.sum(_target_distribution)
    _simulated_distribution = np.array([simulated_distribution["bo4"], 
                                    simulated_distribution["bb"],
                                    simulated_distribution["b5"],
                                    simulated_distribution["b1"],
                                    simulated_distribution["5o4"],
                                    simulated_distribution["ao4"],
                                    simulated_distribution["55"]])
    _simulated_distribution = _simulated_distribution/np.sum(_simulated_distribution)
    
    # Calculate the Euclidean distance
    euclidean_distance = np.linalg.norm(_target_distribution - _simulated_distribution)
    return(euclidean_distance)

def evaluate_minimum_maximum_DP(data):
    """
    Input: data -  A list containing multiple dictionaries. Each dictionary contains information about 1 lignin molecule. 

    Output: The minimum DP observed in a dataset and the maximum DP observed in a dataset
    """

    # Initialize lists to store the values
    values = []

    # Iterate through the list of dictionaries
    for record_dict in data:
        # Check if the 'value' key exists in the dictionary
        if 'DP' in record_dict:
            # Append the 'value' to the values list
            values.append(record_dict['DP'])
    # Find the minimum and maximum values using the min() and max() functions
    if values:
        min_value = min(values)
        max_value = max(values)
        #print("Minimum value:", min_value)
        #print("Maximum value:", max_value)
        return(min_value, max_value)
    else:
        print("No values found for the specified key.")

def exp_func(x, A):
    return np.exp(-A * x)

def line(x, M):
    return(M*x)

def straight_line(x, M, C):
    return(M*x + C)
def exponential_2params(x, A, B):
    y = B*np.exp(-A * x)
    return y

def inverse_exponential(x, A, B):
    y = A*(1-np.exp(-B*x))
    return(y)

def exponential_fit(x, y):
    # Provide an initial guess for parameter A
    initial_guess_A = 0.4
    # Use curve fit with the initial guess for the fitting parameter
    popt, _ = curve_fit(exp_func, x, y, p0=(initial_guess_A))
    A = popt[0]
    return A, popt

def calculate_mean_DP(library):
    """
    Input: The lignin library of one simulated biomass. This will contain structural information about each individual lignin molecule

    Output: Mean value of size distribution
    """
    # Initialize lists to store the values
    values = []
    
    # Iterate through the list of dictionaries
    for _dict in library:
        # Append the 'value' to the values list
        values.append(_dict["DP"])
        #print(values[-1])


    mean_DP = np.mean(values)
    return mean_DP
    

def readin_fitting_results(file_path):

    """
    Input: file_path - The file location of an output file containing each iteration of the fitting parameters in a fitting distribution

    euclidian distance, monomer addition rate, maximum monomer number and minimum monomer number

    Each row contains 1 instance of this. This function reads the final line of the file, indicating this is the end point of the fitting.
    """

    # Open the file in read mode
    with open(file_path, 'r') as file:
        # Read all lines into a list
        lines = file.readlines()

    # Check if the file is not empty
    if lines:
        # Get the last line
        last_line = lines[-1].strip()  # Remove any leading/trailing whitespace

        # Split the last line into four numbers using a delimiter (e.g., space)
        numbers = last_line.split()

        # Check if there are exactly four numbers
        if len(numbers) == 4:
            # Convert the strings to integers or floats as needed
            try:
                euclidian_distance = float(numbers[0])
                monomer_addition_rate = float(numbers[1])
                minimum_monomers = float(numbers[2])
                maximum_monomers = float(numbers[3])

                # Now you have your four numbers
                data = [euclidian_distance, monomer_addition_rate, minimum_monomers, maximum_monomers]
                #print("Data:", data) 
            except ValueError:
                print("Error: Could not convert all elements to numbers.")
        else:
            print("Error: The last line does not contain exactly four numbers.")
        return(data)
    else:
        print("Error: The file is empty.")
    

