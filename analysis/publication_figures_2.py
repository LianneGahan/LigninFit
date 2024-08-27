
import matplotlib.pyplot as plt
import json 
import os 

from rdkit.Chem import AllChem as Chem
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from matplotlib.cm import get_cmap
import pandas as pd
from collections import defaultdict
from rdkit.Chem import Draw
import seaborn as sns
from scipy import stats
from figure_generation_functions import euclidian_dist, calculate_single_fingerprint, compute_similarity_matrix_twolibraries, compute_similarity_matrix, calculate_fingerprints, bond_metrics, line, straight_line, calculate_avg_bond_distribution, calculate_bond_distribution_fitness, exponential_2params, import_library, calculate_mean_DP, evaluate_minimum_maximum_DP, exponential_fit, exp_func, readin_fitting_results
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances_argmin_min
from scipy.stats import mode


# Generate rainbow color palette
cmap = plt.get_cmap('rainbow')

import numpy as np

def load_biomass_data(file_path):
    with open(file_path, 'r') as f:
        data = json.load(f)
    return data

def statistical_tests_bond_distribution(libraries, biomass_labels):
    ####### FOR THE BRANCHING FACTOR
    
    BO4 = np.zeros(len(libraries))
    BB = np.zeros(len(libraries))
    B5 = np.zeros(len(libraries))
    B1 = np.zeros(len(libraries))
    AO4 = np.zeros(len(libraries))
    _4O5 = np.zeros(len(libraries))
    _55 = np.zeros(len(libraries))
    index=0
    labels = ["$\\beta-O4$", "$\\beta-\\beta$", "$\\beta-5$", "$5-5$", "$\\alpha-O4$", "$4-O5$", "$\\beta-1$"]
    for ind, library in libraries:
        bond_distribution = calculate_avg_bond_distribution(library)
        BO4[index] = bond_distribution["bo4"]
        BB[index] = bond_distribution["bb"]
        B5[index] = bond_distribution["b5"]
        B1[index] = bond_distribution["b1"]
        AO4[index] = bond_distribution["ao4"]
        _4O5[index] = bond_distribution["5o4"]
        _55[index] = bond_distribution["55"]
        index += 1
    BONDS = [BO4,BB,B5,_55, AO4, _4O5, B1]
    for bond in BONDS:
        data_bonds = {
            'Group': ['Grasses'] * 4 + ['Hardwoods'] * 7 + ['Softwoods'] * 5 + ["Other"] * 4,
            'bond': bond
        }
        print(np.mean(bond[:4]))
        print(np.mean(bond[4:11]))
        print(np.mean(bond[11:16]))
        print(np.mean(bond[16:]))

        df_bonds = pd.DataFrame(data_bonds)
    # Perform ANOVA
        anova_result = stats.f_oneway(
            df_bonds[df_bonds['Group'] == 'Grasses']['bond'],
            df_bonds[df_bonds['Group'] == 'Hardwoods']['bond'],
            df_bonds[df_bonds['Group'] == 'Softwoods']['bond']
        ) 
        print(f"F-statistic: {anova_result.statistic}, p-value: {anova_result.pvalue}")
        # Perform Tukey's HSD test
        tukey_result = pairwise_tukeyhsd(df_bonds['bond'], df_bonds['Group'], alpha=0.05)
        print("Test results for branching factor")
        print(tukey_result) 



def plot_bond_distributions(libraries, biomass_labels):
    """
    Plots the average bond distribution for each biomass as a stacked bar chart.

    Parameters:
    -----------
    libraries : list
        A list of library objects or data structures containing bond distribution data for different biomasses.
    biomass_labels : list
        A list of strings representing the labels for each biomass, used for the x-axis of the plot.

    Returns:
    --------
    None
        The function saves a stacked bar chart as 'bond_distributions_complete.png', showing the proportion 
        of different bond types for each biomass in the provided libraries.

    Description:
    ------------
    This function creates a stacked bar chart where each bar represents a different biomass and 
    is divided into segments that correspond to the proportions of various bond types, including 
    $\\beta-O4$, $\\beta-\\beta$, $\\beta-5$, $5-5$, $\\alpha-O4$, $4-O5$, and $\\beta-1$. 
    The bond distribution data is calculated for each library using the `calculate_avg_bond_distribution` 
    function. The plot is customized with labels, colors, and font sizes, and is saved as a PNG file.
    """
    fontsize = 22
    fig, ax = plt.subplots(figsize=(16,6))
    # Set the axis labels and title
    ax.set_xlabel('Biomass', fontsize=fontsize)
    ax.set_ylabel('Proportion of bond type', fontsize=fontsize)

    colors = [cmap(i) for i in np.linspace(0, 1, len(target_data)+1)]
    index=0
    # Set positions of each biomass
    x = np.arange(len(libraries))
    ax.set_xticks(x)
    ax.set_xticklabels(biomass_labels)
    transparency = 1
    width = 0.8

    BO4 = np.zeros(len(libraries))
    BB = np.zeros(len(libraries))
    B5 = np.zeros(len(libraries))
    B1 = np.zeros(len(libraries))
    AO4 = np.zeros(len(libraries))
    _4O5 = np.zeros(len(libraries))
    _55 = np.zeros(len(libraries))

    labels = ["$\\beta-O-4$", "$\\beta-\\beta$", "$\\beta-5$", "$5-5$", "$\\alpha-O-4$", "$4-O-5$", "$\\beta-1$"]
    for ind, library in libraries:
        bond_distribution = calculate_avg_bond_distribution(library)
        BO4[index] = bond_distribution["bo4"]
        BB[index] = bond_distribution["bb"]
        B5[index] = bond_distribution["b5"]
        B1[index] = bond_distribution["b1"]
        AO4[index] = bond_distribution["ao4"]
        _4O5[index] = bond_distribution["5o4"]
        _55[index] = bond_distribution["55"]
        index += 1

    colors = [cmap(i) for i in np.linspace(0, 1, 7)]


    ax.bar(x, BO4, width, label=labels[0], edgecolor='white', color = colors[0])
    Bottom = np.zeros(len(BO4))
    for i in range(len(BO4)):
        Bottom[i] = BO4[i]

    ax.bar(x, BB, width, bottom=Bottom, label=labels[1], edgecolor='white', color = colors[1])
    for i in range(len(BO4)):
        Bottom[i] += BB[i]

    ax.bar(x, B5, width, bottom=Bottom, label=labels[2], edgecolor='white', color = colors[2])
    for i in range(len(BO4)):
        Bottom[i] += B5[i]

    ax.bar(x, _55, width, bottom=Bottom, label=labels[3], edgecolor='white', color = colors[3])
    for i in range(len(BO4)):
        Bottom[i] += _55[i]

    ax.bar(x, AO4, width, bottom=Bottom, label=labels[4], edgecolor='white', color = colors[4])
    for i in range(len(BO4)):
        Bottom[i] += AO4[i]
   
    ax.bar(x, _4O5, width, bottom=Bottom, label=labels[5], edgecolor='white', color = colors[5])
    for i in range(len(BO4)):
        Bottom[i] += _4O5[i]
    
    ax.bar(x, B1, width, bottom=Bottom, label=labels[6], edgecolor='white', color = colors[6])

    ax.legend(loc='upper right', ncol=7, bbox_to_anchor=(1.05, -0.17), fontsize=18, markerscale = 1)
    ax.tick_params(axis='x', labelsize=20)  # Change the fontsize of x-axis ticks
    ax.tick_params(axis='y', labelsize=fontsize)  # Change the fontsize of y-axis ticks
    savestring = "bond_distributions_complete.png"
    plt.savefig(savestring,  bbox_inches='tight') 

    return()



def plot_DP_distribution(libraries_list, labels):
    fig, ax = plt.subplots(len(libraries_list), figsize=(12,12), sharex=True)
    fontsize = 22

     # Set the axis labels and title
    ax[-1].set_xlabel('Degree of Polymerisation', fontsize=fontsize)
    ax[2].set_ylabel('Frequency', fontsize=fontsize)
    # Set positions of each biomass
    x_positions = np.arange(len(libraries_list))+1
    colors = [cmap(i) for i in np.linspace(0, 1, len(libraries_list))]

    ax[-1].tick_params(axis='x', labelsize=fontsize)  # Change the fontsize of x-axis ticks
    DPs = []
    transparency = 0.8
    i=0
    for ind, library in libraries_list:
        DP = np.array([data["DP"] for data in library])
        maxi = np.max(DP)
        DPs.append(DP)
        ax[i].tick_params(axis='y', labelsize=fontsize-2)  # Change the fontsize of y-axis ticks

        # Create weighted data by repeating each data point according to its weight
        violin_parts = ax[i].hist(DP, bins=50, label=labels[i], density=True, histtype='step', alpha = transparency, linewidth = 4, color=colors[i])#,showmeans=True)#, colors=colors)
        i+=1

    plt.tight_layout()
    plt.savefig("size_distribution__chosenbiomasses.png", bbox_inches='tight')

def cluster_biomasses(libraries, labels):
    index= 0
    BO4 = np.zeros(len(libraries))
    BB = np.zeros(len(libraries))
    B5 = np.zeros(len(libraries))
    B1 = np.zeros(len(libraries))
    AO4 = np.zeros(len(libraries))
    _4O5 = np.zeros(len(libraries))
    _55 = np.zeros(len(libraries))

    for ind, library in libraries:
        bond_distribution = calculate_avg_bond_distribution(library)
        BO4[index] = bond_distribution["bo4"]
        BB[index] = bond_distribution["bb"]
        B5[index] = bond_distribution["b5"]
        B1[index] = bond_distribution["b1"]
        AO4[index] = bond_distribution["ao4"]
        _4O5[index] = bond_distribution["5o4"]
        _55[index] = bond_distribution["55"]
        index += 1

    data = np.array([BO4, BB, B5, _55, AO4, _4O5, B1])

    # Transpose data to have biomasses as rows
    data = data.T

    # Parameters
    n_clusters = 4
    n_init = 1000  # Number of times to run K-Means
    random_states = np.random.randint(0, 10000, n_init)  # Different random states for each run

    # Store cluster assignments from each run
    all_cluster_assignments = np.zeros((data.shape[0], n_init), dtype=int)

    # Run K-Means multiple times
    for i, random_state in enumerate(random_states):
        kmeans = KMeans(n_clusters=n_clusters, random_state=random_state)
        kmeans.fit(data)
        all_cluster_assignments[:, i] = kmeans.labels_

    # Determine the most common cluster assignment for each biomass
    final_clusters = np.zeros(data.shape[0], dtype=int)
    for i in range(data.shape[0]):
        # Find the most frequent cluster assignment for each biomass
        final_clusters[i] = mode(all_cluster_assignments[i, :])[0][0]

    # Print the resulting cluster assignments
    biomass_clusters = {i: [] for i in range(n_clusters)}
    for biomass_index, cluster_id in enumerate(final_clusters):
        biomass_clusters[cluster_id].append(biomass_index)

    for cluster_id, biomass_indices in biomass_clusters.items():
        print(f"Cluster {cluster_id}: Biomasses {biomass_indices}")
        for i in biomass_indices:
            print(labels[i])

def plot_branching_distribution(libraries_list, labels):
    fig, ax = plt.subplots(figsize=(12,6))
    fontsize = 22

     # Set the axis labels and title
    ax.set_xlabel('Biomass', fontsize=fontsize)
    ax.set_ylabel('DP weighted branching factor', fontsize=fontsize)
    # Set positions of each biomass
    x_positions = np.arange(len(libraries_list))+1
    colors = [cmap(i) for i in np.linspace(0, 1, len(libraries_list))]

    ax.set_xticks(x_positions)
    ax.set_xticklabels(labels, fontsize=fontsize)
    ax.tick_params(axis='x', labelsize=fontsize)  # Change the fontsize of x-axis ticks
    ax.tick_params(axis='y', labelsize=fontsize)  # Change the fontsize of y-axis ticks
    branch = []
    DPs = []
    index =0
    transparency = 1
    for ind, library in libraries_list:
        print(index)
        if index != 30:
            print(labels[index])
            branching = np.array([data["Branches"] for data in library])
            DP = np.array([data["DP"] for data in library])
            DPs.append(DP)
            # Create weighted data by repeating each data point according to its weight
            weighted_data = np.repeat(branching/DP, DP)
            branch.append(weighted_data)
        index+=1
    violin_parts = ax.violinplot(branch,showmeans=True)#, colors=colors)

    # Extract means
    means = [segment[0, 1] for segment in violin_parts['cmeans'].get_segments()]

    # Example data (replace with actual data)
    data = {
        'Group': ['Grasses'] * 4 + ['Hardwoods'] * 6 + ['Softwoods'] * 5,
        'BranchingParameter': means[:15]
    }

    df = pd.DataFrame(data)
# Perform ANOVA
    anova_result = stats.f_oneway(
        df[df['Group'] == 'Grasses']['BranchingParameter'],
        df[df['Group'] == 'Hardwoods']['BranchingParameter'],
        df[df['Group'] == 'Softwoods']['BranchingParameter']
    ) 
    print(f"F-statistic: {anova_result.statistic}, p-value: {anova_result.pvalue}")
    # Perform Tukey's HSD test
    tukey_result = pairwise_tukeyhsd(df['BranchingParameter'], df['Group'], alpha=0.05)
    print(tukey_result)

    # Customize colors
    for i, pc in enumerate(violin_parts['bodies']):
        pc.set_facecolor(colors[i])
        pc.set_edgecolor(colors[i])
        pc.set_alpha(0.7)
    plt.tight_layout()
    plt.savefig("branching_distribution__chosenbiomasses.png", bbox_inches='tight')


def plot_fingerprint_similarity_single_biomass(libraries_list, labels):
    fig, ax = plt.subplots(figsize=(12,12))
    fontsize = 22

     # Set the axis labels and title
    ax.set_xlabel('Biomass', fontsize=fontsize)
    ax.set_ylabel('DP weighted branching factor', fontsize=fontsize)
    for ind, lib in libraries_list:
        library = lib
    fingerprints = calculate_fingerprints(library)
    similarity_matrix = compute_similarity_matrix(fingerprints)
    print(np.shape(similarity_matrix))
        # Check if similarity_matrix is empty or contains only NaN values
    if similarity_matrix.size == 0:# or np.isnan(similarity_matrix).all():
        raise ValueError("The similarity matrix is empty or contains only NaN values")

    #similarity_matrix = pairwise_distances(fingerprints, metric='jaccard')
    sns.heatmap(similarity_matrix, cmap='viridis')
    plt.title('Similarity Matrix Heatmap')
    plt.savefig("fingerprintsim.png", bbox_inches='tight')

def plot_fingerprint_similarity(libraries_list, biomasses, labels):
    fig, ax = plt.subplots(figsize=(15,12))
    fontsize = 22
    index = 0
    biomasses_fingerprints = []
    # Sort libs by sg ratio



    print("Sorting library data and calculating fingerprints")
    for ind, lib in libraries_list:
        print(index)
        fingerprints = []

        data = [data for data in lib if data["DP"]>=10 and data["DP"]<=10]
        euclidians = [bond_metrics(biomasses[index], mol["Bonds"])["euclidian"] for mol in data]
        # Combine data with their corresponding euclidean distances
        data_with_euclidians = list(zip(data, euclidians))

        # Sort the combined list based on euclidean distances
        sorted_data_with_euclidians = sorted(data_with_euclidians, key=lambda x: x[1])

        # Extract the sorted data from the sorted list of tuples
        sorted_data = [item[0] for item in sorted_data_with_euclidians]

        for i in range(1):
            try:
                fingerprint = calculate_single_fingerprint(sorted_data[i])
                print(sorted_data_with_euclidians[i][1])
                fingerprints.append(fingerprint)

            except:
                fingerprint = [0]
        biomasses_fingerprints.append(fingerprints)
        print("There are :", len(fingerprints))
        index+=1
    n = len(biomasses)
    similarity_matrix = np.zeros((n, n))
    print("Calculating similarities")
    for i in range(n):
        for j in range(i, n):
            similarity_matrix_pair = compute_similarity_matrix_twolibraries(np.array(biomasses_fingerprints[i]), np.array(biomasses_fingerprints[j]))
            similarity_matrix[i, j] = np.mean(similarity_matrix_pair)
            similarity_matrix[j, i] = np.mean(similarity_matrix_pair)
        # Check if similarity_matrix is empty or contains only NaN values
    if similarity_matrix.size == 0:# or np.isnan(similarity_matrix).all():
        raise ValueError("The similarity matrix is empty or contains only NaN values")

    """ax.axhline(4, color="white", linestyle="--", linewidth=2)
    ax.axhline(11, color="white", linestyle="--", linewidth=2)
    ax.axhline(16, color="white", linestyle="--", linewidth=2)
    ax.axvline(4, color="white", linestyle="--", linewidth=2)
    ax.axvline(11, color="white", linestyle="--", linewidth=2)
    ax.axvline(16, color="white", linestyle="--", linewidth=2)"""

    ticks = [0.5, 0.75, 1.0]
    #similarity_matrix = pairwise_distances(fingerprints, metric='jaccard')
    sns.heatmap(similarity_matrix, cmap='Greens', cbar_kws={'ticks': ticks, "label": "Tanimoto similarity"})
     # Set the axis labels and title    
    cbar = ax.collections[0].colorbar

    cbar.set_label('Tanimoto similarity', fontsize=fontsize)

    ax.collections[0].colorbar.set_ticklabels([abs(t) for t in ticks], fontsize=fontsize)
    plt.gca().invert_yaxis()

    x_positions = np.arange(len(libraries_list))+0.5

    ax.set_xlabel('Biomass', fontsize=fontsize)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(labels, fontsize=fontsize-1)
    ax.set_yticks(x_positions)
    ax.set_yticklabels(labels, fontsize=fontsize-1)#, rotation = 90)
 
    ax.set_ylabel('Biomass', fontsize=fontsize)
    #plt.title('Similarity Matrix Heatmap')
    plt.savefig("biomass_fingerprints_comparison_sgratioorder.png", bbox_inches='tight')

def plot_average_bond_distribution_similarity(libraries, biomass_labels):
    """
    Generates a heatmap of Spearman correlation coefficients comparing the bond distributions across different biomass libraries.

    Parameters:
    -----------
    libraries : list
        A list of tuples, where each tuple contains metadata (ignored in this function) and a list of molecules.
        Each molecule is represented as a dictionary with bond type counts.
    biomass_labels : list of str
        A list of labels corresponding to each biomass library, used for axis labels in the heatmap.

    Functionality:
    --------------
    This function:
    - Calculates the average bond distribution for each biomass library using a helper function `calculate_avg_bond_distribution`.
    - Computes the Spearman correlation coefficient between the bond distributions of each pair of biomass libraries.
    - Creates a heatmap to visualize the correlation matrix, with customized colors, tick labels, and axis labels.
    - Saves the heatmap as 'bond_distribution_similarity.png'.

    Plot Customization:
    -------------------
    - The heatmap uses the 'Blues_r' colormap and displays tick labels for specific correlation values.
    - Font sizes and transparency are set for better visualization.
    - Axis labels and tick labels are set according to the provided biomass labels.

    Returns:
    --------
    The heatmap is saved as a PNG file ('bond_distribution_similarity.png') in the current directory.
    """

    fontsize = 22
    fig, ax = plt.subplots(figsize=(12,10))
    # Set the axis labels and title


    colors = [cmap(i) for i in np.linspace(0, 1, len(target_data)+1)]
    index=0
    # Set positions of each biomass
    x = np.arange(len(libraries))
    ax.set_xticks(x)
    ax.set_xticklabels(biomass_labels)
    transparency = 1
    width = 0.8
    outputs = []
    for i in range(len(libraries)):
        bond_distribution = calculate_avg_bond_distribution(libraries[i][1])
        bonds = bond_distribution.values()
            #print(bonds)
        outputs.append(list(bonds)[:])


    outputs = np.array(outputs)

    columns = biomass_labels#["$\\beta$-O-4", "$\\beta-\\beta$", "$\\beta$-5", "5-5", "$\\alpha$-O-4", "4-O-5", "$\\beta$-1"]
    #columns = ['A', 'B', 'C', 'W', 'X', 'Y', 'Z', 'A', 'B', 'C', 'W', 'X', 'Y', 'Z']
    #import pandas as pd
    #df = pd.DataFrame(outputs, columns=columns)
    correlation_matrix = np.zeros([len(libraries), len(libraries)])
    for i in range(len(libraries)):
        for j in range(len(libraries)):
            _distribution1 = outputs[i, :]
            print(_distribution1)
        
            _distribution2 = outputs[j, :]
        
            # Normalise so we keep the ratios of bond types
            _distribution1 = _distribution1/np.sum(_distribution1)
            _distribution2 = _distribution2/np.sum(_distribution2)
            correlation_matrix[i,j] = stats.spearmanr(_distribution1, _distribution2).statistic


    #correlation_matrix = df.corr(method='pearson')

    ticks = [0.5, 0.6, 0.7, 0.8, 0.9, 1]

    sns.heatmap(correlation_matrix, cmap='Blues_r', cbar_kws={'ticks': ticks, "label": "Spearman Correlation"})#, vmin=-1, vmax=1)
    cbar = ax.collections[0].colorbar
    ax.collections[0].colorbar.set_ticklabels([t for t in ticks], fontsize=fontsize)
    cbar.set_label('Spearman Correlation', fontsize=fontsize)
    ax.set_xlabel('Biomass', fontsize=fontsize+2)
    ax.set_ylabel('Biomass', fontsize=fontsize+2)
    ax.set_xticklabels(biomass_labels)
    ax.set_yticklabels(biomass_labels)#, rotation=90)


    # Compute the correlation matrix
    ax.tick_params(axis='x', labelsize=fontsize-2)  # Change the fontsize of x-axis ticks
    ax.tick_params(axis='y', labelsize=fontsize-2)#, rotation=90)  # Change the fontsize of y-axis ticks
    plt.gca().invert_yaxis()
    plt.tight_layout()
    savestring = "bond_distribution_similarity.png"
    plt.savefig(savestring, bbox_inches='tight')


def plot_moleculespread(libraries_list, target_data, labels):
    fig, ax = plt.subplots(figsize=(12, 6))
    fontsize = 22

    # Set the axis labels and title
    ax.set_xlabel('Biomass', fontsize=fontsize)
    ax.set_ylabel('Euclidean distance', fontsize=fontsize)
    # Set positions of each biomass
    x_positions = np.arange(1, len(libraries_list)+1)

    ax.set_xticks(x_positions)
    ax.set_xticklabels(labels, fontsize=fontsize)
    ax.tick_params(axis='x', labelsize=fontsize)  # Change the fontsize of x-axis ticks
    ax.tick_params(axis='y', labelsize=fontsize)  # Change the fontsize of y-axis ticks
    euclidians = []
    averages = []

    index = 0
    averages =[]
    for ind, library in libraries_list:
        bond_distribution= calculate_avg_bond_distribution(library=library)
        print(bond_distribution, target_data[index],)
        averages.append(euclidian_dist(target_data[index], bond_distribution))
        ed = []
        print(target_data[index])
        for data in library:
            # Check if at least one value is non-zero
            at_least_one_non_zero = any(value != 0 for value in data["Bonds"].values())
            if at_least_one_non_zero and data["DP"] > 1:
                ed.append(euclidian_dist(target_data[index], data["Bonds"]))
        ed = np.array(ed)

        # Flatten the array and remove NaNs and infs
        ed_flat = ed.flatten()
        ed_cleaned = ed_flat[~np.isnan(ed_flat) & ~np.isinf(ed_flat)]

        euclidians.append(ed_cleaned)
        index += 1

    violin_parts = ax.violinplot(euclidians, showextrema=True)

    # Customize colors and lines
    for i, pc in enumerate(violin_parts['bodies']):
        pc.set_facecolor("grey")
        pc.set_edgecolor("grey")
        pc.set_alpha(0.7)


    # Customize range lines (extrema) if they exist
    if 'cmins' in violin_parts:
        violin_parts['cmins'].set_color('black')
    if 'cmaxes' in violin_parts:
        violin_parts['cmaxes'].set_color('black')
    if 'cbars' in violin_parts:
        violin_parts['cbars'].set_color('black')



    # Add average lines inside each violin plot
    print(averages)
    for i in range(len(averages)):
            print(i)
            x = np.array([x_positions[i] - 0.4, x_positions[i] + 0.4])  # Adjust for the width of each violin
            y = np.array([averages[i], averages[i]])
            ax.plot(x, y, color='k', linestyle='--', linewidth=3)#, label=f'Average Euclidean Distance (Biomass {labels[i]})')

    # Add legend
    handles, legend_labels = ax.get_legend_handles_labels()
    unique_labels = list(dict.fromkeys(legend_labels))  # Remove duplicate labels
    unique_handles = [handles[legend_labels.index(label)] for label in unique_labels]

    plt.tight_layout()
    plt.savefig("euclidiandistance_spread.png", bbox_inches='tight')


def plot_molecule(smiles, name):
    """
    Generates a 2D depiction of a molecule from a SMILES string and saves it as an SVG and PNG file.

    Parameters:
    -----------
    smiles : str
        The SMILES (Simplified Molecular Input Line Entry System) string representing the molecule to be drawn.
    name : str
        The desired name for the output PNG file, including the extension (e.g., 'molecule.png').

    Functionality:
    --------------
    This function:
    - Converts the SMILES string into a RDKit molecule object.
    - Creates a 2D SVG drawing of the molecule with customized visual options:
    - Larger bond line width for clarity.
    - Stereo annotation for chiral centers.
    - Increased font size for atom labels.
    - Black and white color palette for a classic look.
    - Saves the SVG drawing to 'output.svg'.
    - Converts the SVG file to a PNG file using the `cairosvg` library for easier viewing.

    Requirements:
    -------------
    - RDKit: For generating the molecule and drawing it.
    - cairosvg: For converting the SVG file to a PNG format.

    Returns:
    --------
    None

    The molecule is saved as both an SVG file ('output.svg') and a PNG file with the specified name.
    """

    # Generate molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    # Assuming 'mol' is your molecule and 'name' is your desired output file name
    # Assuming 'mol' is your molecule and 'name' is your desired output file name
    d2d = Draw.MolDraw2DSVG(1800, 600)  # SVG with smaller size to check scaling
    d2d.drawOptions().bondLineWidth = 2.0  # Thicker bonds for clarity
    d2d.drawOptions().addStereoAnnotation = True
    d2d.drawOptions().atomLabelFontSize = 100  # Larger font size for atom labels
    d2d.drawOptions().minFontSize = 22
    d2d.drawOptions().useBWAtomPalette()  # Use black and white palette for a classic look

    # Draw the molecule
    d2d.DrawMolecule(mol)
    d2d.FinishDrawing()

    # Save the SVG image
    # Save the SVG image
    svg = d2d.GetDrawingText()
    with open("output.svg", "w") as f:
        f.write(svg)

    # Optionally, convert SVG to PNG for easier viewing
    # This requires additional libraries like cairosvg
    import cairosvg
    cairosvg.svg2png(url="output.svg", write_to=name)

    print("Molecule saved")

def plot_highlightmolecules():
    smiles = "[CH3][O][c]1[cH][c]([CH]2[O][c]3[c]([O][CH3])[cH][c]([CH]([OH])[CH]([CH2][OH])[O][c]4[c]([O][CH3])[cH][c]([CH]([OH])[CH]([CH2][OH])[O][c]5[c]([O][CH3])[cH][c]([CH]([OH])[CH]([CH2][OH])[O][c]6[c]([O][CH3])[cH][c]([CH]7[O][CH2][CH]8[CH]([c]9[cH][c]([O][CH3])[c]([OH])[c]([O][CH3])[cH]9)[O][CH2][CH]78)[cH][c]6[O][CH3])[cH][c]5[O][c]5[c]([O][CH3])[cH][c]([CH]6[O][CH2][CH]7[CH]([c]8[cH][c]([O][CH3])[c]([O][CH]([CH2][OH])[CH]([OH])[c]9[cH][c]([O][CH3])[c]([O][CH]([CH2][OH])[CH]([OH])[c]%10[cH][c]([O][CH3])[c]([O][CH]([CH2][OH])[CH]([OH])[c]%11[cH][c]([O][CH3])[c]([O][CH]([CH2][OH])[CH]([OH])[c]%12[cH][c]([O][CH3])[c]([O][CH]([CH2][OH])[CH]([OH])[c]%13[cH][cH][c]([O][CH]([CH2][OH])[CH]([OH])[c]%14[cH][cH][c]([O][CH]([CH2][OH])[CH]([OH])[c]%15[cH][c]([O][CH3])[c]([O][CH]([CH2][OH])[CH]([OH])[c]%16[cH][c]([O][CH3])[c]([O][CH]([CH2][OH])[CH]([OH])[c]%17[cH][c]([O][CH3])[c]([O][CH]([CH2][OH])[CH]([OH])[c]%18[cH][c]([O][CH3])[c]([O][CH]([CH2][OH])[CH]([OH])[c]%19[cH][cH][c]([O][CH]([CH2][OH])[CH]([OH])[c]%20[cH][c]([O][CH3])[c]%21[c]([cH]%20)[CH]([CH2][OH])[CH]([c]%20[cH][c]([O][CH3])[c]([OH])[c]([O][CH3])[cH]%20)[O]%21)[c]([O][CH3])[cH]%19)[c]([O][CH3])[cH]%18)[c]([O][CH3])[cH]%17)[c]([O][CH3])[cH]%16)[c]([O][CH3])[cH]%15)[c]([O][CH3])[cH]%14)[c]([O][CH3])[cH]%13)[c]([O][CH3])[cH]%12)[c]([O][CH3])[cH]%11)[c]([O][CH3])[cH]%10)[c]([O][c]%10[c]([O][CH3])[cH][c]([CH]([OH])[CH]([CH2][OH])[O][c]%11[c]([O][CH3])[cH][c]([CH]%12[O][CH2][CH]%13[CH]([c]%14[cH][c]([O][CH3])[c]([OH])[c]([O][CH3])[cH]%14)[O][CH2][CH]%12%13)[cH][c]%11[O][CH3])[cH][c]%10[O][CH3])[cH]9)[c]([O][CH3])[cH]8)[O][CH2][CH]67)[cH][c]5[O][CH3])[cH][c]4[O][CH3])[cH][c]3[CH]2[CH2][OH])[cH][cH][c]1[OH]"

    # Generate molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)

    # Add hydrogens to the molecule
    mol_with_h = Chem.AddHs(mol)

    # Define SMARTS patterns for guaiacyl and syringyl units
    guaiacyl_smarts = "COc1ccccc1O"  # Guaiacyl unit: -OCH3 group on a phenol ring
    syringyl_smarts = "COc1cc(OC)c(O)cc1"  # Syringyl unit: two -OCH3 groups on a phenol ring

    # Convert SMARTS patterns to molecule objects
    guaiacyl_pattern = Chem.MolFromSmarts(guaiacyl_smarts)
    syringyl_pattern = Chem.MolFromSmarts(syringyl_smarts)

    # Find all substructure matches
    guaiacyl_matches = mol_with_h.GetSubstructMatches(guaiacyl_pattern)
    syringyl_matches = mol_with_h.GetSubstructMatches(syringyl_pattern)

    # Get the atom indices to highlight
    highlight_atoms_guaiacyl = set(atom for match in guaiacyl_matches for atom in match)
    highlight_atoms_syringyl = set(atom for match in syringyl_matches for atom in match)

    # Set up drawing options
    d2d = Draw.MolDraw2DCairo(1000, 1000)
    d2d.drawOptions().bondLineWidth = 2.0
    d2d.drawOptions().addStereoAnnotation = True
    d2d.drawOptions().atomLabelFontSize = 14
    d2d.drawOptions().useBWAtomPalette()
    d2d.drawOptions().includeAtomNumbers = False  # Optionally disable atom numbering

    # Define a function to remove hydrogen atoms from the molecule
    def remove_hydrogens(mol):
        """Remove hydrogen atoms from the molecule."""
        return Chem.RemoveHs(mol)

    # Remove hydrogen atoms from the molecule
    mol_without_h = remove_hydrogens(mol_with_h)

    # Draw the molecule with highlighted atoms
    highlight_colors = {}
    for atom in highlight_atoms_guaiacyl:
        highlight_colors[atom] = (1, 0, 0)  # Red color for guaiacyl units
    for atom in highlight_atoms_syringyl:
        highlight_colors[atom] = (0, 0, 1)  # Blue color for syringyl units

    d2d.DrawMolecule(mol_without_h, highlightAtoms=highlight_colors.keys(), highlightAtomColors=highlight_colors)
    d2d.FinishDrawing()


    # Save the drawing to a PDF file
    d2d.WriteDrawingText('molecule_guaiacyl_syringyl_highlighted.png')

def plot_DP_bonddistributions(library, label):

    """
    Plots the distribution of bond types as a function of degree of polymerization (DP) for a given library of molecules.

    Parameters:
    -----------
    library : tuple
        A tuple where the first element is typically ignored and the second element is a list of dictionaries.
        Each dictionary represents a molecule and contains:
        - "DP": int, degree of polymerization of the molecule.
        - "Bonds": dict, a dictionary where keys are bond types (e.g., 'bo4', 'bb', etc.) and values are their counts.
    
    label : str
        A label used for the plot title or file naming, though not currently implemented in the function.

    Functionality:
    --------------
    This function:
    - Groups molecules into bins based on their degree of polymerization (DP) in increments of 2.
    - Sums the occurrences of different bond types within each DP bin.
    - Normalizes these sums to get the proportion of each bond type within each bin.
    - Plots a stacked bar chart representing the distribution of bond types across different DP bins.

    The x-axis represents the bin index, and the y-axis represents the proportion of each bond type within each bin.
    Bond types are displayed in a specific order ('bo4', 'bb', 'b5', '55', 'ao4', '5o4', 'b1') and are colored using
    a colormap ('rainbow').

    The plot is saved as "DP_bonddist.png" in the current directory.

    Notes:
    ------
    - The function adjusts for bins that might have no data by assigning them a proportion of zero.
    - The `label` parameter is included in the function signature but is not currently used within the function.

    Returns:
    --------
    None
    """
    fig, ax = plt.subplots(figsize=(16, 6))
    fontsize = 24

    ax.set_xlabel('Degree of polymerisation', fontsize=fontsize)
    ax.set_ylabel('Proportion of bond type', fontsize=fontsize)
    ax.set_ylim(0, 1.05)
    ax.tick_params(axis='x', labelsize=fontsize, length=5, width=2)
    ax.tick_params(axis='y', labelsize=fontsize)

    molecules = library[1]
    
    bin_size = 2
    dp_bins = [(i, i + bin_size - 1) for i in range(1, 111, bin_size)]
    
    binned_bond_sums = defaultdict(lambda: defaultdict(int))

    for molecule in molecules:
        dp = molecule["DP"]
        bonds = molecule["Bonds"]
        for (bin_min, bin_max) in dp_bins:
            if bin_min <= dp <= bin_max:
                for bond_type, count in bonds.items():
                    binned_bond_sums[(bin_min, bin_max)][bond_type] += count

    normalized_distributions = {}
    for dp_bin, bond_sums in binned_bond_sums.items():
        total_bonds = sum(bond_sums.values())
        if total_bonds > 0:
            normalized_distributions[dp_bin] = {bond: count / total_bonds for bond, count in bond_sums.items()}
        else:
            normalized_distributions[dp_bin] = {bond: 0 for bond in bond_sums.keys()}

    bond_order = ['bo4', 'bb', 'b5', '55', "ao4", '5o4', 'b1']
    bond_types = [bond for bond in bond_order if any(bond in dist for dist in normalized_distributions.values())]

    cmap = get_cmap('rainbow')
    colors = [cmap(i) for i in np.linspace(0, 1, 7)]

    latex_labels = {"bo4":"$\\beta-O-4$","bb": "$\\beta-\\beta$","b5": "$\\beta-5$","55": "$5-5$", "ao4": "$\\alpha-O-4$", "5o4": "$4-O-5$", "b1": "$\\beta-1$"}

    color_mapping = {bond: colors[i] for i, bond in enumerate(bond_types)}
    
    bottoms = np.zeros(len(dp_bins))

    for bond_type in bond_types:
        x_values = list(range(len(dp_bins)))
        y_values = [normalized_distributions.get(dp_bin, {}).get(bond_type, 0) for dp_bin in dp_bins]
        
        ax.bar(x_values, y_values, bottom=bottoms, edgecolor='white', color=color_mapping[bond_type], label=latex_labels[bond_type])
        bottoms += y_values

    # Set custom x-axis labels
    dp_labels = [bin_max for bin_min, bin_max in dp_bins]
    tick_positions = [0] + list(np.arange(4, len(dp_bins), 5))  # Include the first position (0)
    tick_labels = ['1'] + [str(dp_labels[i]) for i in tick_positions[1:]]
    
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels)

    ax.legend(loc='upper right', ncol=7, fontsize=fontsize-8, bbox_to_anchor=(0.9, -0.25))
    plt.tight_layout()

    plt.savefig("DP_bonddist.png", bbox_inches="tight")

def plot_mass_comparison(libraries_list, target_data, labels):
    fig, ax = plt.subplots(figsize=(12,6))
    fontsize = 22


     # Set the axis labels and title
    ax.set_xlabel('Biomass', fontsize=fontsize)
    ax.set_ylabel('DP', fontsize=fontsize)
    # Set positions of each biomass
    x_positions = np.arange(len(libraries_list))+1
    colors = [cmap(i) for i in np.linspace(0, 1, len(libraries_list))]

    ax.set_xticks(x_positions)
    ax.set_xticklabels(labels, fontsize=fontsize)
    ax.tick_params(axis='x', labelsize=fontsize)  # Change the fontsize of x-axis ticks
    ax.tick_params(axis='y', labelsize=fontsize)  # Change the fontsize of y-axis ticks
    masses = []
    DPs=[]
    index =0
    transparency=1
    S_monolignol_weight =210
    G_monolignol_weight = 180
    for ind, library in libraries_list:
        DP = np.array([data["DP"] for data in library])
        DPs.append(DP)
        DP = np.mean(DP)#calculate_mean_DP(library)
        if 'MW' in target_data[index]:
            target_DP = (target_data[index]["MW"]*(target_data[index]["Monolignols"]["S"]/100)/S_monolignol_weight) + (target_data[index]["MW"]*(target_data[index]["Monolignols"]["G"]/100)/G_monolignol_weight)
            if index==0:
                ax.scatter(x_positions[index], target_DP, marker="v", alpha=transparency, label = "Experimental Data", color = "k", s=300)
                #ax.scatter(x_positions[index], DP, marker="v", alpha=transparency,  label = "Simulation Data",  color = "k", s=300)
            else:
                ax.scatter(x_positions[index], target_DP, marker="v", alpha=transparency,  color = "k", s=300)
            #ax.scatter(x_positions[index], DP, marker="v", alpha=transparency,  color = "k", s=300)
        index+=1
    violin_parts = ax.violinplot(DPs,showmeans=True)#, colors=colors)
    # Customize colors
    for i, pc in enumerate(violin_parts['bodies']):
        pc.set_facecolor(colors[i])
        pc.set_edgecolor(colors[i])
        pc.set_alpha(0.7)
    fig.legend(fontsize=18, loc='lower center', bbox_to_anchor=(0.5, -0.12), markerscale = 1, ncols=2)

    plt.tight_layout()
    plt.savefig("DPcomparison__chosenbiomasses.png", bbox_inches='tight')


def order_data_by_cluster():
    # Desired order of labels in clusters
    desired_order = ['H2', 'H3', 'H6', 'H7', 'S1', 'X1', 'X2', 'X3', 'X4',  # Cluster 0
                     'G1', 'G2', 'G3', 'G4', 'H5', 'S2',                     # Cluster 1
                    'H4', 'S3', 'S4', 'S5',                                  # Cluster 2
                    'H1']                                                    # Cluster 3

    # Load biomass data
    biomasses = load_biomass_data("biomasses.json")

    # Create a mapping from biomass label to its position in the desired order
    label_order_map = {label: i for i, label in enumerate(desired_order)}

    # Sort the biomass labels based on the desired order
    sorted_dict_names = sorted(biomasses.keys(), key=lambda x: label_order_map[biomasses[x]["short label"]])

    location_string = "example_data/set_C/"

    target_data = []
    libs = []
    fitting_data = []

    for biomass_label in sorted_dict_names:
        library_import_string = location_string + biomass_label + "/BEST_FIT/best_Run/Output/library.json"
        if os.path.exists(library_import_string):
            libs.append(import_library(library_import_string, "ligninkmc_"))
            fitting_import_string = location_string + biomass_label + "/best_kin_specs.txt"
            fitting_data.append(np.loadtxt(fitting_import_string))
            target_data.append(biomasses[biomass_label])
    return(libs, target_data)

biomasses = load_biomass_data("biomasses.json")
def get_sort_key(dict_name):
    return biomasses[dict_name]["short label"]

def get_sg_key(dict_name):
    return biomasses[dict_name]["sg_ratio"]
sorted_dict_names = sorted(biomasses.keys(), key=get_sort_key)

location_string = "example_data/set_C/"

target_data = []

allparamslimitedsglibs =[]

fitting_params_list = []
fitting_data = []
for biomass_label in sorted_dict_names:
    library_import_string = location_string + biomass_label + "/BEST_FIT/best_Run/Output/library.json"
    if os.path.exists(library_import_string):
        allparamslimitedsglibs.append(import_library(library_import_string, "ligninkmc_"))
        fitting_import_string = location_string + biomass_label + "/best_kin_specs.txt"
        fitting_data.append(np.loadtxt(fitting_import_string))
        target_data.append(biomasses[biomass_label])

    else: print(library_import_string)
fitting_params_list.append(fitting_data)




labels = [dict["short label"] for dict in target_data]
#plot_branching_distribution(allparamslimitedsglibs, labels)
"""]plot_DP_distribution(allparamslimitedsglibs, labels)
plot_mass_comparison(allparamslimitedsglibs, target_data, labels)"""
plot_moleculespread(allparamslimitedsglibs, target_data, labels)




#plot_average_bond_distribution_similarity(allparamslimitedsglibs, labels)


#### FOR FINGERPRINTS
sorted_dict_names = sorted(biomasses.keys(), key=get_sg_key)
location_string = "example_data/set_C/"
target_data = []
allparamslimitedsglibs =[]
fitting_params_list = []
fitting_data = []

for biomass_label in sorted_dict_names:
    library_import_string = location_string + biomass_label + "/BEST_FIT/best_Run/Output/library.json"
    if os.path.exists(library_import_string):
        allparamslimitedsglibs.append(import_library(library_import_string, "ligninkmc_"))
        fitting_import_string = location_string + biomass_label + "/best_kin_specs.txt"
        fitting_data.append(np.loadtxt(fitting_import_string))
        target_data.append(biomasses[biomass_label])
labels = [dict["short label"] for dict in target_data]
#plot_fingerprint_similarity(allparamslimitedsglibs, target_data,labels)
#cluster_biomasses(allparamslimitedsglibs, labels)
#plot_DP_bonddistributions(allparamslimitedsglibs[0], "G1")
#statistical_tests_bond_distribution(allparamslimitedsglibs, labels)
#plot_highlightmolecules()
def process_lib_for_plotmolecule():
    poplar_lib = allparamslimitedsglibs[8][1]
    DPs = np.array([data["DP"] for data in poplar_lib])
    poplar_smile = np.array([data["smilestring"] for data in poplar_lib])
    for i, smile in enumerate(poplar_smile):
        if DPs[i] < 25 and DPs[i] > 15:
            name = "birchmols/birch_mol__" + str(i) + ".png"
            plot_molecule(poplar_smile[i], name)
    return()

"""smile = np.array([data["smilestring"] for data in allparamslimitedsglibs[8][1]])
dp = np.array([data["DP"] for data in allparamslimitedsglibs[8][1]])
sgs = np.array([data["sg_ratio"] for data in allparamslimitedsglibs[8][1]])
bonds = np.array([data["Bonds"] for data in allparamslimitedsglibs[8][1]])


plot_molecule(smile[422], "birch_molecule.png")"""

"""
libs, td = order_data_by_cluster()
labs = [dict["short label"] for dict in td]

print(labs)
plot_bond_distributions(libs, labs)"""