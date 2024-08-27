from __future__ import print_function
from common_wrangler.common import (MAIN_SEC, GOOD_RET, INPUT_ERROR, KB, PLANCK_CONST_JS, KCAL_MOL_TO_J_PART,
                                    INVALID_DATA, OUT_DIR, InvalidDataError, warning, process_cfg, make_dir,
                                    create_out_fname, str_to_file, round_sig_figs)
from ligninkmc.kmc_functions import (run_kmc, generate_mol, find_fragments, fragment_size, break_bond_type, gen_tcl)
from ligninkmc.create_lignin import (calc_rates, create_initial_monomers, create_initial_events,
                                     create_initial_state, analyze_adj_matrix,  adj_analysis_to_stdout)
from ligninkmc.kmc_common import (DEF_E_BARRIER_KCAL_MOL, ADJ_MATRIX, MONO_LIST, MONOMER, OX, GROW, C, Monomer, Event)
from ligninkmc.kmc_common import (BO4, B5, BB, B1, B1_ALT, C5O4, AO4, C5C5, G, S, C,
                                   ADJ_MATRIX, BONDS, CHAIN_LEN, RCF_YIELDS)

from scipy.sparse import dok_matrix
from scipy.optimize import curve_fit
import cProfile
import pstats

from rdkit.Chem.Draw import MolToFile
from rdkit.Chem.rdMolInterchange import MolToJSON
from rdkit.Chem import AllChem, Descriptors
from rdkit import Chem

import matplotlib.pyplot as plt
import json
import numpy as np
#Parallelization
import multiprocessing

# For performance
import time


temp = 298.15  # K
rxn_rates = calc_rates(temp, ea_kcal_mol_dict=DEF_E_BARRIER_KCAL_MOL)

TCL_NAME = "psfgen.tcl"
PSF_NAME = 'lignin'
TOPPAR_DIR = "toppar/"

#fun = par.delayed(run_kmc)


def gaussian(x, amplitude, mean, stddev):
    return amplitude * np.exp(-((x - mean) / stddev) ** 2 / 2)

def fit_gaussian_histogram(data, histogram, bin_edges):
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    initial_guess = [1.0, np.mean(data), np.std(data)]
    params, covariance = curve_fit(gaussian, bin_centers, histogram, p0=initial_guess)

    amplitude_fit, mean_fit, stddev_fit = params
    print(mean_fit, amplitude_fit, stddev_fit)
    fitted_curve = gaussian(bin_centers, amplitude_fit, mean_fit, stddev_fit)
    
    #filtered_indices = np.where((data >= lower_bound) & (data <= upper_bound))[0]
    #filtered_data = data[filtered_indices]
    return(bin_centers, fitted_curve, mean_fit, stddev_fit, amplitude_fit)
def calculate_bond_distribution(data_dictionary):
    """
    Calculate the normalised distribution of specific bond types within a collection of simulated structures.

    Inputs:
        simulation_result (dictionary): A collection of simulation data containing adjacency matrices and bond information.
        num_sims (int): The number of simulations or structures in the provided data.

    Returns:
        dict: A dictionary containing the normalized distribution of bond types within the simulated structures.
            The keys are bond type names, and the values are the normalized counts of each bond type."""
    #get all bond dictionaries
    Bonds = [inner_dict["Bonds"] for inner_dict in data_dictionary.values() if "Bonds" in inner_dict]
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
        #print(dictionary)
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


def reassign_unique_ids(list_of_list_of_dicts):
    """
    Reassigns unique identifiers to a list of dictionaries and compiles them into a single dictionary.

    This function takes a list of lists, where each inner list contains dictionaries. It assigns
    a unique ID to each dictionary in the format "ligninkmc_X", where X is an incrementing integer, 
    and compiles all the dictionaries into a single output dictionary.

    Parameters:
    -----------
    list_of_list_of_dicts : list of lists of dict
        A list containing sublists, each with dictionaries representing data entries.

    Returns:
    --------
    new_dictionary : dict
        A dictionary where each key is a unique ID in the format "ligninkmc_X" and each value 
        is a dictionary from the input list.
    """

    unique_id = 1  # Start with ID 1
    new_dictionary = {}
    for s in list_of_list_of_dicts:

        for num, dict in enumerate(s):
            #print(s[num])
            new_dictionary["ligninkmc_" + str(unique_id)] = s[num]
            unique_id += 1  # Increment unique ID
    return new_dictionary

def split_molecule_with_aromaticity(molecule):
    """
    Splits a molecule into fragments while ensuring correct aromaticity perception.

    This function takes an RDKit molecule object, adjusts aromaticity perception, and then
    splits the molecule into separate fragments.

    Parameters:
    -----------
    molecule : RDKit.Chem.Mol
        An RDKit molecule object representing the chemical structure.

    Returns:
    --------
    molecules : tuple of RDKit.Chem.Mol
        A tuple containing RDKit molecule objects, each representing a fragment of the original molecule.
    """
    # Ensure aromaticity is perceived correctly
    Chem.SanitizeMol(molecule, sanitizeOps=Chem.SANITIZE_ALL^Chem.SANITIZE_SETAROMATICITY)

    # Split the molecule into separate molecules (fragments)
    molecules = Chem.GetMolFrags(molecule, asMols=True)
    return molecules


def SMILES_pipeline(adjacency_matrix, monomer_list):
    """
    Generates SMILES strings for each molecule in a lignin polymer from an adjacency matrix and a list of monomers.

    Parameters:
    adjacency_matrix (dict of tuple of int): A dictionary where the keys are tuples representing pairs of monomers, 
                                             and the values represent the bond types between these monomers. For example, 
                                             `((i, j): bond_type)` where `bond_type` can be BO4, B1, BB, etc.
    monomer_list (list of Monomer): A list where each element is a `Monomer` object. Each `Monomer` object should have
                                     attributes necessary to construct the molecular structure (e.g., type).

    Returns:
    dict: A dictionary where each key is an integer index corresponding to a molecule. Each value is another dictionary 
          containing the following properties for that molecule:
          - "smilestring": The SMILES string representation of the molecule.

    Description:
    This function performs the following steps to generate SMILES strings for lignin molecules:
    1. Uses the `generate_mol` function to create an RDKit molecule block from the provided adjacency matrix and monomer list.
    2. Converts the RDKit molecule block to an RDKit molecule object.
    3. Splits the molecule into individual fragments or components using aromaticity considerations.
    4. Converts each fragment into a SMILES string using RDKit's SMILES writing functionality.
    5. Stores the SMILES strings in a dictionary with indices as keys and the SMILES strings as values.

    Note:
    - Ensure that the `generate_mol` function correctly generates a molecular block in a format compatible with RDKit.
    - The `split_molecule_with_aromaticity` function should correctly handle aromatic and non-aromatic structures.
    - RDKit's `Chem` module functions (`MolFromMolBlock`, `MolToSmiles`, etc.) require RDKit to be installed and properly configured in your environment.
    """
    # Default out is SMILES, which requires getting an rdKit molecule object; also required for everything
    #    except the TCL format
    dict={}
    block = generate_mol(adjacency_matrix, monomer_list)
    molecule_total = Chem.MolFromMolBlock(block)
    molecules = split_molecule_with_aromaticity(molecule_total)
    #Save SMILES string for individual molecule
    params = Chem.SmilesWriteParams()
    params.allHsExplicit = True
    params.canonical = True
    params.doIsomericSmiles = True
    smiles = [Chem.MolToSmiles(molecule, params) for molecule in molecules]
    #Write dictionaryoutput to return to workflow
    for i in range(len(molecules)):
        dict[i] = {
                "smilestring": smiles[i],}
                #"MW": mw}
    return(dict)

def generate_analysis_parallel(adjacency_matrix, monomer_list, file_writing, savefile_path):
    """
    Analyzes the bond distribution and generates SMILES strings for a given lignin molecule structure.

    This function takes an adjacency matrix representing the lignin molecule, a list of monomers,
    and performs a complete workflow to analyze the bond distribution and calculate the degree of polymerization (DP).
    It then generates SMILES strings for the molecule and associates them with the respective bond information.

    Parameters:
    -----------
    adjacency_matrix : np.ndarray
        A square matrix representing the connectivity between atoms in the lignin molecule, where
        each entry indicates the presence or absence of a bond.
    
    monomer_list : list
        A list of monomers in the lignin molecule that correspond to the indices in the adjacency matrix.
    
    file_writing : bool
        A flag indicating whether to write the output data to a file.
    
    savefile_path : str
        The path where the output file should be saved if file_writing is set to True.
    
    Returns:
    --------
    lignin_dictionaries : list of dict
        A list of dictionaries where each dictionary contains bond information and the corresponding
        SMILES string for a segment of the lignin molecule.

    Notes:
    ------
    - The function uses the `separate_molecule_properties` function to calculate bond information
      for each part of the lignin molecule.
    - The `SMILES_pipeline` function generates the corresponding SMILES strings for each segment of the molecule.
    - The function handles cases where B-1 bonds are present by concatenating the SMILES strings of bonded segments.
    - If `file_writing` is True, the function saves the results to the specified file path.
    """
    lignin_dictionaries = []

    # Complete work flow to calculate bond distribution and DP for each lignin molecule
    partial_dict_A = separate_molecule_properties(adjacency_matrix, monomer_list)
    #print("Bond Analysis complete!")

    partial_dict_B = SMILES_pipeline(adjacency_matrix, monomer_list)
    
    smiles_dictionary_id = 0
    bond_dictionary_id = 0
    for i in range(len(partial_dict_A)):
    # Cycling through bond info dictionary, but need to check smiles dict if B-1 bonds are present
        if partial_dict_A[bond_dictionary_id]['Bonds']["b1"] != 0:    
            smiles = ".".join([partial_dict_B[smiles_dictionary_id + k]["smilestring"] for k in range(partial_dict_A[bond_dictionary_id]['Bonds']["b1"])])
            smiles_dictionary_id += partial_dict_A[bond_dictionary_id]['Bonds']["b1"]
        else:
            smiles = partial_dict_B[smiles_dictionary_id]["smilestring"]
        
        # Combine dictionaries
        partial_dict_A[bond_dictionary_id]["smilestring"] = smiles
        bond_dictionary_id+=1
        smiles_dictionary_id+=1
    #print("Analysis complete!")

    #combined dictionaries
    return(partial_dict_A)

def save_structures(id_list, molecules, dict):
    """
    Saves individual molecule images of type ".png" and ".svg".
    Saves a SMILES string of the molecule
    """

    return()

def separate_molecule_properties(adjacency_matrix, monomer_list):
    """
    Analyzes molecular properties from an adjacency matrix and a list of monomers, and returns a dictionary containing
    detailed information about each lignin molecule.

    Parameters:
    adjacency_matrix (dict of tuple of int): A dictionary where the keys are tuples representing pairs of monomers, and
                                             the values represent the bond types between these monomers. For example, 
                                             `((i, j): bond_type)` where `bond_type` can be BO4, B1, BB, etc.
    monomer_list (list of Monomer): A list where each element is a `Monomer` object. Each `Monomer` object should have
                                     a `type` attribute that specifies the type of the monomer (e.g., 'syringyl', 'guaiacol').

    Returns:
    dict: A dictionary where each key represents a unique molecule identified by its index. Each value is another dictionary
          containing the following properties for that molecule:
          - "Bonds": A dictionary with bond types as keys and their counts as values.
          - "sg_ratio": The S/G ratio of the molecule, calculated as the ratio of syringyl to guaiacol monomers.
          - "DP": The degree of polymerization, which is the number of monomers in the molecule.
          - "Branches": The number of branches in the molecule.
          - "Monolignols": A dictionary with counts of different types of monolignols ('S', 'G', 'C') in the molecule.

    Description:
    This function processes an adjacency matrix and a list of monomers to extract and analyze individual lignin molecules.
    It performs the following steps:
    1. Uses the adjacency matrix to find fragments (molecules) and the number of branches in each molecule.
    2. For each molecule, initializes a bond count dictionary to keep track of different bond types.
    3. Constructs an adjacency matrix specific to the molecule and calculates the number of each bond type based on the
       global adjacency matrix.
    4. Computes the S/G ratio of the molecule and identifies the number of branches.
    5. Compiles the information into a dictionary that provides detailed properties for each molecule.

    Note:
    - Ensure that the `find_fragments` function is defined and correctly returns the connected components and branch counts.
    - The `Monomer` class or type used in `monomer_list` should have a `type` attribute to identify the monomer type.
    - The bond types in `bonding_dict` (BO4, B1, BB, etc.) should be defined and consistent with those used in `adjacency_matrix`.
    """

    bonding_dict = {(4, 8): BO4, (8, 4): BO4, (8, 1): B1, (1, 8): B1, (8, 8): BB, (5, 5): C5C5,
                    (8, 5): B5, (5, 8): B5, (7, 4): AO4, (4, 7): AO4, (5, 4): C5O4, (4, 5): C5O4}
    dict = {}


    #Provides tuples full of monomer indexes, each tuple contains monomers for one lignin molecule. Also provides number of branches in each molecule
    connected, branches = find_fragments(adjacency_matrix)
    # for each molecule
    #Looping over each molecule 
    
    #print(len(connected))
    for number, molecule in enumerate(connected):
        #print("Molecule Analysis complete!")
        molecule = list(molecule)


        #initialise the dictionary  
        bond_count_dict = {BO4: 0,  BB: 0, B5: 0, B1: 0, C5O4: 0, AO4: 0, C5C5: 0}

        #make an empty dok matrix to copy values over to from the global adj matrix
        molecule_adj_matrix = dok_matrix((len(molecule), len(molecule)), dtype=int)
        #make an empty monomer list to copy values over from the global monolist
        #to calculate SG ratio properties
        s_lignol_sum, g_lignol_sum, c_lignol_sum = 0, 0, 0

        for i in range(0, len(molecule)):
            # Confirm the identity of the monomer    
            if monomer_list[molecule[i]].type == 'syringyl':
                s_lignol_sum += 1
            elif monomer_list[molecule[i]].type == 'guaiacol':
                g_lignol_sum += 1
            else: c_lignol_sum += 1

            for j in range(i+1, len(molecule)):
                #print(i, j)
                mol1 = molecule[i]
                mol2 = molecule[j]
                #copy value of adjacency matrix at index mol1 and mol2
                molecule_adj_matrix[i, j] = adjacency_matrix[(mol1, mol2)]
                #Check number at location of ij in adj matrix
                if i != j:
                    #Access bond type by looking at bound carbon location in adjacency matrix
                    bond = adjacency_matrix[(mol1, mol2)], adjacency_matrix[(mol2, mol1)]

                    #so long as the values are not zero, match them according to the bonding dictionary
                    if all(bond):
                        #add one count to the correct dictionary
                        bond_count_dict[bonding_dict[bond]] += 1


        #calculate sg ratio of molecule:
        if g_lignol_sum == 0:
            sg_ratio = s_lignol_sum
        elif s_lignol_sum == 0:
            sg_ratio = 0
        else: sg_ratio = s_lignol_sum/g_lignol_sum

        #Identify number of branches that a particular lignin molecule has
        branching_number = branches[number]

        dict[number] = {
            "Bonds": bond_count_dict,
            "sg_ratio": sg_ratio,
            "DP": len(molecule),
            "Branches": branching_number,
            "Monolignols": {"S": s_lignol_sum,
                            "G": g_lignol_sum,
                            "C": c_lignol_sum},
        }
    return dict

def adjust_energy_barriers(energy_barrier_dict, scale_factor_GG, scale_factor_SS, scale_factor_SG):
    """
    Adjusts the energy barriers in a dictionary based on provided scaling factors for different types of monomer pairs.

    Parameters:
    energy_barrier_dict (dict): A dictionary where the keys represent different reaction types and the values are nested dictionaries.
                                The nested dictionaries contain energy barrier values for different monomer pairs.
    scale_factor_GG (float): The scaling factor to be applied to energy barriers associated with the G-G monomer pair.
    scale_factor_SS (float): The scaling factor to be applied to energy barriers associated with the S-S monomer pair.
    scale_factor_SG (float): The scaling factor to be applied to energy barriers associated with the S-G or G-S monomer pairs.

    Returns:
    dict: The updated `energy_barrier_dict` with adjusted energy barrier values based on the provided scaling factors.

    Description:
    This function iterates over the `energy_barrier_dict` to adjust the energy barriers for different monomer pairs. 
    The dictionary is expected to have a structure where each key maps to a nested dictionary that contains energy barriers 
    for different monomer pairs. The function applies the following scaling factors:
    - `scale_factor_GG` to energy barriers for the G-G monomer pair.
    - `scale_factor_SS` to energy barriers for the S-S monomer pair.
    - `scale_factor_SG` to energy barriers for both S-G and G-S monomer pairs.

    Note:
    - Ensure that the keys in the nested dictionaries match the expected monomer pair tuples (e.g., (G, G), (S, G)).
    - The function assumes that energy barrier values are numeric (integers or floats) and will be multiplied by the corresponding scale factors.

    Example:
    adjusted_energies = adjust_energy_barriers(energy_barriers, scale_factor_GG=1.1, scale_factor_SS=0.9, scale_factor_SG=1.05)

    """
    # Iterate over each key in the dictionary
    for key, value in energy_barrier_dict.items():
        # Check if the value is a dictionary
        if isinstance(value, dict):
            # Iterate over each key in the nested dictionary
            for nested_key, nested_value in value.items():
                if nested_key == (G, G):  
                    for subnested_key, subnested_value in nested_value.items():
                        # Assuming subnested_value is an integer or float
                        value[nested_key][subnested_key] *= scale_factor_GG
                if nested_key == (S, G): 
                    for subnested_key, subnested_value in nested_value.items():
                        # Assuming subnested_value is an integer or float
                        value[nested_key][subnested_key] *= scale_factor_SG
                if nested_key == (G, S): 
                    for subnested_key, subnested_value in nested_value.items():
                        # Assuming subnested_value is an integer or float
                        value[nested_key][subnested_key] *= scale_factor_SG
                if nested_key == (S, S): 
                    for subnested_key, subnested_value in nested_value.items():
                        # Assuming subnested_value is an integer or float
                        value[nested_key][subnested_key] *= scale_factor_SS
    return(energy_barrier_dict)

def generate_analysis(results, num_sims, num_cores, file_writing, savefile_path):
    """
    Analyzes the results of a series of lignin-KMCK simulations using parallel processing.
    This function calculates bond distributions, SMILES structures, branching coefficients, and chemical functions
    from the simulation outputs.

    Parameters:
    results (list of dict): The raw results from the KMC simulations, where each dictionary in the list contains 
                            information such as adjacency matrices and monomer lists for each simulation.
    num_sims (int): The number of simulation repeats used in the analysis.
    num_cores (int): The number of CPU cores to be used for parallel processing.
    file_writing (bool): Flag indicating whether to write the analysis results to a file.
    savefile_path (str): The file path where the JSON results will be saved if `file_writing` is True.

    Returns:
    dict: A dictionary containing the final analyzed data, which includes information on bond distributions, 
          SMILES structures, branching coefficients, and chemical functions.

    Description:
    This function processes the raw output from KMC simulations by performing the following steps:
    1. Extracts adjacency matrices and monomer lists from the simulation results.
    2. Utilizes multiprocessing to parallelize the analysis across multiple CPU cores. 
    3. Calls the `generate_analysis_parallel` function to compute the required metrics for each simulation.
    4. Flattens the resulting dictionaries and reassigns unique IDs using `reassign_unique_ids`.
    5. Optionally writes the final results to a JSON file if `file_writing` is set to True.

    Notes:
    - Ensure that `generate_analysis_parallel` and `reassign_unique_ids` functions are correctly implemented and available in the context where this function is called.
    - Make sure to handle exceptions or errors related to file operations if `file_writing` is True.
    - The `multiprocessing.Pool` object requires that all functions used within the pool be pickleable.

    Example:
    ```python
    final_data = generate_analysis(results, num_sims=10, num_cores=4, file_writing=True, savefile_path="output.json")
    ```
    """

    #Simulation needs to be split into different adj matrices for different structures 
    #Also need monomer lists from each simulation
    #adjacency matrix for each simulation
    cur_adjs = [results[j][ADJ_MATRIX] for j in range(num_sims)] 
    #list of monomers in each simulation
    monolist = [results[j][MONO_LIST] for j in range(num_sims)]
    #assign an id to each simulation 
    #simulations = ["ligninkmc_" + str(i) for i in range(len(results))]

    # Multiprocessing code is directly placed here
    with multiprocessing.Pool(processes=num_cores) as pool:
        # Lignin libraries is a list of list of dictionaries
        lignin_libraries = pool.starmap(generate_analysis_parallel, [(cur_adjs[i], monolist[i], file_writing, savefile_path) for i, result in enumerate(results)])

    # Let's flatten the dictionary and reassign the IDs
    final_lignin_dictionary = reassign_unique_ids(lignin_libraries)

    #if file_writing:
    print("Dictionaries built!") 
    print("Converting and writing to json file...")
    json_data = json.dumps(final_lignin_dictionary)
    with open(savefile_path, "w") as outfile:
        outfile.write(json_data)

    print("Your work here is done!")
    return(final_lignin_dictionary)


def simulation(sg_ratio, ini_num_monos, max_num_monos, mono_add_rate, simulation_reaction_rates, t_max):
    """
    Runs a Kinetic Monte Carlo (KMC) simulation for lignin polymerization.

    Parameters:
    sg_ratio (float): The ratio of S to G monolignols in the simulation. This ratio is used to determine the initial distribution of S and G monomers.
    ini_num_monos (int): The initial number of monolignols to be generated for the simulation.
    max_num_monos (int): The maximum number of monolignols allowed in the simulation.
    mono_add_rate (float): The rate at which new monolignols are added to the system.
    simulation_reaction_rates (dict): A dictionary containing the reaction rates for various processes in the simulation.
    t_max (float): The maximum simulation time.

    Returns:
    result: The outcome of the KMC simulation, which may include information such as the final state of the system, reaction counts, or time evolution data.

    Description:
    This function sets up and runs a Kinetic Monte Carlo (KMC) simulation for lignin polymerization based on the provided parameters. It begins by determining the percentage of S monolignols based on the provided sg_ratio and generates an initial set of monomers. It then initializes the simulation state, sets up initial events, and runs the KMC simulation until the maximum number of monomers or the maximum simulation time is reached.

    The simulation output typically includes the final state of the system, reaction statistics, or other relevant metrics as defined by the KMC implementation.

    Note:
    - Ensure that the `create_initial_monomers`, `create_initial_events`, `create_initial_state`, and `run_kmc` functions are defined and implemented correctly for this simulation function to work.
    - The `Event` class or function used in `initial_events.append` should be properly defined and imported.
    """
    pct_s = sg_ratio / (1 + sg_ratio)
    # Make choices about what kinds of monomers there are and create them
    # Use a random number and the given sg_ratio to determine the monolignol types to be initially modeled
    monomer_draw = np.random.rand(ini_num_monos)
    initial_monomers = create_initial_monomers(pct_s, monomer_draw)

    # Initialize the monomers, events, and state
    initial_events = create_initial_events(initial_monomers, simulation_reaction_rates)
    initial_state = create_initial_state(initial_events, initial_monomers)
    initial_events.append(Event(GROW, [], rate=mono_add_rate))
    result = run_kmc(simulation_reaction_rates, initial_state, initial_events, n_max = max_num_monos, t_max = t_max, 
                                            sg_ratio = sg_ratio)
    return(result)



