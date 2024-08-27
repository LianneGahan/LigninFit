
import matplotlib.pyplot as plt
import json 
from rdkit.Chem import AllChem as Chem
from rdkit import DataStructs
from sklearn.metrics import r2_score
from adjustText import adjust_text
from scipy.optimize import curve_fit
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from scipy import stats
import matplotlib.gridspec as gridspec
from collections import Counter, OrderedDict
import pandas as pd
from statsmodels.stats.multicomp import pairwise_tukeyhsd

from rdkit.Chem import Draw, MolFromSmiles, Descriptors
from scipy.sparse import dok_matrix
import matplotlib.patches as patches
import seaborn as sns

#from rdkit.Chem import Fragments
import os 
import re

from figure_generation_functions import calculate_fingerprints, bond_metrics, line, straight_line, calculate_avg_bond_distribution, calculate_bond_distribution_fitness, exponential_2params, import_library, calculate_mean_DP, evaluate_minimum_maximum_DP, exponential_fit, exp_func, readin_fitting_results
from matplotlib import colors
# Generate rainbow color palette
cmap = plt.get_cmap('rainbow')

import numpy as np




def plot_compare_SGratio_freeenergies__heatmap(fitting_parameters_list, biomasses, biomass_labels):
    """
    This function creates a multi-panel plot to visualize changes in specific fitting parameters across various 
    biomass samples. The plot includes four subplots:
    1. Logarithm of S-S scaling parameter vs biomass.
    2. Logarithm of G-G scaling parameter vs biomass.
    3. Logarithm of S-G scaling parameter vs biomass.
    4. Logarithm of monomer addition rate vs biomass.

    Parameters:
    -----------
    fitting_parameters_list : list
        A list of fitting parameter sets, where each set corresponds to a different fitting scenario. Each fitting 
        parameter set contains the scaling parameters for different types of bonds and the monomer addition rate.
    
    biomasses : list of dicts
        A list containing information about the biomasses, including their types (e.g., Softwood, Hardwood, etc.).
    
    biomass_labels : list of str
        A list of labels corresponding to each biomass sample.

    Returns:
    --------
    None
        The function generates and saves a multi-panel plot comparing the fitting parameters for different biomasses.

    Notes:
    ------
    The function automatically adjusts layout and positioning of elements to ensure clarity. A legend is included 
    at the bottom center of the plot to differentiate between different fitting scenarios.
    """

    fontsize = 22
    fig = plt.figure(figsize=(16,12))

# Create a figure

    # Create a GridSpec object with 5 rows and 1 column, where the 4th row will be a spacer
    gs = gridspec.GridSpec(4, 1, height_ratios=[1, 1, 1, 1])# 0.0, 1])  # Adjust height_ratios as needed
    # Add subplots to the figure
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)
    ax3 = fig.add_subplot(gs[2])
    ax4 = fig.add_subplot(gs[3])

    ax = [ax1, ax2, ax3, ax4]
     # Set the axis labels and title
    #ax[1].set_ylabel('Scaling parameter change from Original', fontsize=fontsize)
    cmap = plt.get_cmap('gnuplot2')
    colors = [cmap(i) for i in np.linspace(0, 1, 12)]
    color1 = colors[4]
    color2 = colors[7]
    color3 = colors[8]


    # Set positions of each biomass

    ax[0].tick_params(axis='x', labelsize=fontsize)  # Change the fontsize of x-axis ticks
    ax[1].tick_params(axis='x', labelsize=fontsize)  # Change the fontsize of x-axis ticks
    ax[2].tick_params(axis='x', labelsize=fontsize)  # Change the fontsize of x-axis ticks
    ax[0].tick_params(axis='y', labelsize=fontsize)  # Change the fontsize of y-axis ticks
    ax[1].tick_params(axis='y', labelsize=fontsize)  # Change the fontsize of y-axis ticks
    ax[2].tick_params(axis='y', labelsize=fontsize)  # Change the fontsize of y-axis ticks
    ax[3].tick_params(axis='x', labelsize=fontsize)  # Change the fontsize of x-axis ticks
    ax[3].tick_params(axis='y', labelsize=fontsize)  # Change the fontsize of y-axis ticks

    transparency = 0.8
    scalings = []
    for i, fitting_parameters in enumerate(fitting_parameters_list):
        gg_scaler = [np.log(data[4]) for data in fitting_parameters]
        ss_scaler = [np.log(data[5]) for data in fitting_parameters]
        sg_scaler = [np.log(data[6]) for data in fitting_parameters]
        sg_ratios = [data[3] for data in fitting_parameters]
        monaddrates = [data[0] for data in fitting_parameters]

        scalings.append([gg_scaler, ss_scaler, sg_scaler, monaddrates, sg_ratios])

    labs = ["ln(S-S)",  "ln(G-G)", "ln(S-G)", "log(Monomer\naddition rate)"]
    x_positions = np.arange(len(scalings[0][4][:]))

        #    gg_scaler = [data[4] for data in fitting_parameters]
    colors = [color3, color2,  color1, ]
    barindex =0
    ax[0].axhline(0, color=color1, linestyle="--", linewidth=2)
    ax[1].axhline(0, color=color1, linestyle="--", linewidth=2)
    ax[2].axhline(0, color=color1, linestyle="--", linewidth=2)    
    ax[3].set_ylabel(labs[3], fontsize=fontsize)
    #ax[3].
    icons = ["v", "s", "o"]
    #for j in range(len(scalings[0][0][:])):
    for i in range(3):
        ax[i].set_ylabel(labs[i], fontsize=fontsize)

        for j in range(len(scalings[i][4][:])):
            if biomasses[j]["Type"] == "Softwood":
                icon = "v"
            if biomasses[j]["Type"] == "Hardwood":
                icon = "v"
            if biomasses[j]["Type"] == "Other":
                icon = "v"
            if biomasses[j]["Type"] == "Grass":
                icon = "v"

            
        #label = labs[i]#
            #ax[0].scatter(scalings[i][4][j], scalings[i][0][j], marker=icon, alpha=transparency, color = colors[i], s=300)
            #ax[1].scatter(scalings[i][4][j], scalings[i][1][j], marker=icon, alpha=transparency, color = colors[i], s=300)
            #ax[2].scatter(scalings[i][4][j], scalings[i][2][j], marker=icon, alpha=transparency, color = colors[i], s=300)
            ax[1].scatter(j, scalings[i][0][j], marker=icons[i], alpha=transparency, linewidth=4, facecolor='none', edgecolors=colors[i], s=300)
            ax[0].scatter(j, scalings[i][1][j], marker=icons[i], alpha=transparency,  linewidth=4,  facecolor='none',edgecolors=colors[i], s=300)
            ax[2].scatter(j, scalings[i][2][j], marker=icons[i], alpha=transparency, linewidth=4,   facecolor='none',edgecolors=colors[i], s=300)
            ax[3].scatter(j, scalings[i][3][j], marker=icons[i], alpha=transparency, linewidth=4,  facecolor='none',edgecolors=colors[i], s=300)

    ax[0].scatter([], [], marker="o",  facecolor='none', edgecolor = color1,  linewidth=4, zorder=1, s=300, label = "A: Monomer addtion rate and S/G ratio free")

    ax[0].scatter([], [], marker="s",  facecolor='none',edgecolor = color2,  linewidth=4, zorder=1, s=300, label = "B: All parameters free")

    ax[0].scatter([], [], marker="v",  facecolor='none', edgecolor = color3,  linewidth=4,zorder=1,  s=300, label = "C: All parameters free, limited S/G")

    # Set ticks for each subplot
    x_positions = np.arange(len(scalings[0][4][:]))
    for axes in ax:
        axes.set_xticks(x_positions)
        
    # Hide x-tick labels for the first three subplots
    for i in range(3):
        ax[i].set_xticklabels([])


    ax[1].set_ylim(-6, 3)
    ax[0].set_ylim(-3.5, 3.5)
    ax[2].set_ylim(-4.5, 1.5)
    ax[3].set_ylim(-0.5, 5)


    # Share x-axis for the first three subplots
    for axes in ax:#[:2]:
        axes.get_shared_x_axes().join(axes, ax[0])
        axes.set_xticklabels([])
    #ax[2].set_xticks(x_positions)
    #ax[2].set_xticklabels(biomass_labels)
    #ax[2].set_xlabel('Biomass Labels', fontsize=fontsize)
    ax[3].set_xticks(x_positions)
    ax[3].set_xticklabels(biomass_labels)
    ax[3].set_xlabel('Biomass Labels', fontsize=fontsize)
    #fig.plot([],[], label = "G-G Interactions")
    fig.legend( loc='lower center',bbox_to_anchor=(0.5, -0.07),ncols=2, fontsize=fontsize) 

    # Adjust the layout to prevent overlapping titles
    #fig.tight_layout()

    savestring = "energybarrier_change_sgratios.png"
    fig.savefig(savestring,  bbox_inches='tight') 



def load_biomass_data(file_path):
    with open(file_path, 'r') as f:
        data = json.load(f)
    return data


def plot_compare_SG_ratios_barchart(library_list, biomasses, labels):
    """
    Plots a comparative bar chart of S/Gratios for different experimental conditions against the expected S/Gratio from biomass data.

    Parameters:
    -----------
    library_list : list of lists
        A list containing sublists of dictionaries where each sublist corresponds to a different simulation condition. 
        Each dictionary in the sublist contains bond distribution data and S/Gratio data for various biomasses.
        The structure is assumed to be: `library_list[condition][biomass_index][data_dict]`, where `data_dict` has keys "sg_ratio" and "DP".
        
    biomasses : list of dicts
        A list of dictionaries where each dictionary contains the experimental SG ratio data for each biomass. 
        Each dictionary must have a key "sg_ratio" to access the SG ratio value.

    labels : list of str
        A list of strings representing the labels for each biomass. These labels will be used for the x-axis tick labels.

    Returns:
    --------
    None
        This function does not return any value. It creates and saves a bar chart as an image file ('compare_fit_sgratios_barchart.png').
        
    Notes:
    ------
    - The function compares SG ratios across three different simulation conditions (labelled A, B, C) and plots these against 
      the corresponding SG ratio from experimental biomass data.
    - The bar chart displays the SG ratios for each biomass and condition, with an additional bar for the experimental SG ratio.
    - The function assumes that the `calculate_avg_bond_distribution` and `bond_metrics` functions are defined elsewhere and used 
      to compute average bond distributions and bond metrics, respectively.
    - The colors and labels for the bars are hardcoded, with colors being chosen from a 'gnuplot2' colormap.
    - The plot includes custom text labels below each set of bars and displays a legend for clarity.
    """

    fontsize = 22
    fig, axes = plt.subplots(figsize=(15,4))
     # Set the axis labels and title
    axes.set_xlabel('Biomass', fontsize=fontsize)
    axes.set_ylabel('S/G ratio', fontsize=fontsize)
    # Bar position and shape properties 
    dist = 0.0
    subtr = -0.05
    separation = 0.2
    width = 0.65  # Reduce the width for the target distribution
    barwidth = 0.22
    cmap = plt.get_cmap('gnuplot2')
    colors = [cmap(i) for i in np.linspace(0, 1, 12)]
    color1 = colors[4]
    color2 = colors[7]
    color3 = colors[8]
    barindex=0
    index = 0
    color4 = colors[2]
    colorlist= [color1, color2, color3]
    # Set positions of each biomass
    x_positions = np.arange(len(library_list[0]))
    axes.set_xticks(x_positions-1 + 0.25)
    axes.set_xticklabels(labels)

    plotlablist = ["A", "B", "C"]

    #new_euclidian_distances =  np.zeros(len(new_libraries))
    #old_euclidian_distances =  np.zeros(len(old_libraries))
    transparency = 1
    for i in range(len(library_list[0])):
        bond_distribution_ALLFREE = calculate_avg_bond_distribution(library_list[0][index][1])
        bond_distribution_ALLFREELIMSG = calculate_avg_bond_distribution(library_list[1][index][1])
        bond_distribution_MONADDSG = calculate_avg_bond_distribution(library_list[2][index][1])
        bond_dist_list = [bond_distribution_MONADDSG, bond_distribution_ALLFREE, bond_distribution_ALLFREELIMSG, ]
        
        
        sgratios = [data["sg_ratio"] for data in library_list[0][index][1]]
        DPs = [data["DP"] for data in library_list[0][index][1]]
        sg_ALLFREE = np.average(sgratios, weights = DPs)

        # Obtain all of the SG ratios of the different molecules
        sgratios = [data["sg_ratio"] for data in library_list[1][index][1]]
        DPs = [data["DP"] for data in library_list[1][index][1]]
        sg_ALLFREELIMSG = np.average(sgratios, weights = DPs)
            # Obtain all of the SG ratios of the different molecules
        sgratios = [data["sg_ratio"] for data in library_list[2][index][1]]
        DPs = [data["DP"] for data in library_list[2][index][1]]
        sg_MONADDSG = np.average(sgratios, weights = DPs)

        sg_list = [sg_MONADDSG, sg_ALLFREE, sg_ALLFREELIMSG, ]
        fitness = np.zeros(3)
        for sim in range(len(library_list)):
            bond_distribution = bond_dist_list[sim]
            fitness[sim] = round(bond_metrics(target_data[index], bond_distribution)["euclidian"], 2)
            #axes.text(barindex-1+ separation*sim + barwidth + dist, sg_list[sim] + 0.6, str(fitness[sim]), rotation=90, fontsize=11,horizontalalignment = 'center')
            axes.text(barindex-1+ separation*sim + barwidth + dist, -0.8, plotlablist[sim], fontsize=12,horizontalalignment = 'center')
            axes.bar(barindex-1+ separation*sim + barwidth + dist, sg_list[sim], barwidth,  edgecolor='white', color=colorlist[sim])
        expe_sg =biomasses[index]["sg_ratio"]
        if expe_sg <= 0.1:
            expe_sg = 0.1
        axes.bar(barindex-1 , expe_sg, barwidth - subtr,  edgecolor='white', color="k")
        index +=1
        barindex +=1


    # Define the legend labels and colors
    legend_labels = ["Experimental data","A: Monomer addition rate and S/G ratio free", "B: All parameters free", "C: All parameters free, limited S/G"]
    legend_colors = ["k", color1, color2, color3]

    # Create custom legend handles using Patch for stacked bars
    legend_handles = [patches.Patch(color=color, label=label) for color, label in zip(legend_colors, legend_labels)]



    #bond_distribution_MONADDSG, bond_distribution_ALLFREE, bond_distribution_ALLFREELIMSG,
    # Create a legend for the main plot
    legend = fig.legend(handles=legend_handles, labels=legend_labels, loc='lower center', bbox_to_anchor=(0.5, -0.25), fontsize=fontsize, ncol=2)
    savestring = "compare_fit_sgratios_barchart.png"
    # Adjust the layout to prevent overlapping titles    
    axes.tick_params(axis='x', labelsize=fontsize)  # Change the fontsize of x-axis ticks
    axes.tick_params(axis='y', labelsize=fontsize)  # Change the fontsize of y-axis ticks
    axes.set_ylim(-1, 12.5)
    plt.tight_layout()
    plt.savefig(savestring, bbox_inches='tight')


def plot_compare_bond_distributions(library_list, target_data, biomass_labels):
    """
    Plots the average bond distribution and experimentally determined target distribution
    for each of the biomasses as a stacked bar chart. The plot consists of three subplots:
    one for the complete bond distribution, one for a partial distribution including more bond types,
    and one for a minimal partial distribution.

    Parameters:
    -----------
    library_list : list of lists
        A list containing sublists, each with bond distribution data for different biomass samples.
        Each sublist contains bond distributions from different simulation conditions.
    
    target_data : list of dicts
        A list where each dictionary contains the experimentally determined target bond distribution 
        for each biomass sample, along with metadata such as the level of bond detail.
        Keys include:
        - 'bo4': Proportion of β-O-4 bonds.
        - 'bb': Proportion of β-β bonds.
        - 'b5': Proportion of β-5 bonds.
        - 'b1': Proportion of β-1 bonds.
        - 'ao4': Proportion of α-O-4 bonds.
        - '4o5': Proportion of 4-O-5 bonds.
        - '55': Proportion of 5-5 bonds.
        - 'bond_detail_level': Level of bond detail ('Complete', 'terrell', 'klose').
    
    biomass_labels : list of str
        A list of labels corresponding to the different biomass samples, used for labeling the x-axis.

    Returns:
    --------
    None
        This function creates and saves a plot as 'bond_distributions_and_targets.png', 
        but does not return any values.

    Notes:
    ------
    - The function creates three subplots:
        1. Complete experimental bond distribution.
        2. Partial distribution with more bond types ($\\beta-O-4$, $\\beta-\\beta$, $\\beta-5$, $\\beta-1$, and grouped).
        3. Minimal partial distribution ($\\beta-O-4$, $\\beta-\\beta$, and $\\beta-5$).
    - It also adds a custom legend explaining the colors used for each bond type.
    - The bar heights in each subplot represent the proportions of different bond types for both 
      the target (experimental) and simulated bond distributions.
    """
    fontsize = 20
    fig, ax = plt.subplots(3, figsize=(16, 11), sharey=True)


    # Set the axis labels and title
    ax[1].set_ylabel('Proportion of Bond Type', fontsize=fontsize)
    ax[2].set_xlabel('Biomass', fontsize=fontsize)

    ax[0].set_title("Complete experimental bond distribution", fontsize=fontsize, loc="center")

    ax[1].set_title("Partial experimental bond distribution ($\\beta-O-4$, $\\beta-\\beta$, $\\beta-5$, $\\beta-1$ and grouped)", fontsize=fontsize, loc="center")

    ax[2].set_title("Partial experimental bond distribution ($\\beta-O-4$, $\\beta-\\beta$ and $\\beta-5$)", fontsize=fontsize, loc="center")
    #ax[2].set_ylabel('Proportion of Bond Type', fontsize=fontsize)
    cmap = plt.get_cmap('rainbow')

    colors = [cmap(i) for i in np.linspace(0, 1, 7)]
    index = 0
    # Bar position and shape properties 
    dist = 0.08
    subtr = -0.05
    separation = 0.2
    width = 0.65  # Reduce the width for the target distribution
    barwidth = 0.17
    # To store labels for each subplot
    ax0labels = []
    ax1labels = []
    ax2labels = []
    # To determine axis positions on each subplot
    ax0index = 0 
    ax1index = 0
    ax2index = 0
    separation = 0.2
    width = 0.65  # Reduce the width for the target distribution
    barwidth = 0.17


    plotlablist = ["A", "B", "C"]
    labels = ["$\\beta-O-4$", "$\\beta-\\beta$", "$\\beta-5$", "$5-5$", "$\\alpha-O-4$", "$4-O-5$", "$\\beta-1$"]
    for i in range(len(library_list[0])):
        target_BO4 = target_data[index]["bo4"]
        target_BB = target_data[index]["bb"]
        target_B5 = target_data[index]["b5"]
        target_B1 = target_data[index]["b1"]
        target_AO4 = target_data[index]["ao4"]
        target_4O5 = target_data[index]["4o5"]
        target_55 = target_data[index]["55"]
        target_OTHER = target_data[index]["ao4"] + target_data[index]["4o5"] +target_data[index]["55"]

        BO4 = np.zeros(3)
        BB = np.zeros(3)
        B5 =np.zeros(3)
        B1 = np.zeros(3)
        AO4 = np.zeros(3)
        _4O5= np.zeros(3)
        _55 = np.zeros(3)
        OTHER =np.zeros(3)

        bond_distribution_ALLFREE = calculate_avg_bond_distribution(library_list[0][index][1])
        bond_distribution_ALLFREELIMSG = calculate_avg_bond_distribution(library_list[1][index][1])
        bond_distribution_MONADDSG = calculate_avg_bond_distribution(library_list[2][index][1])

        bond_dist_list = [bond_distribution_MONADDSG, bond_distribution_ALLFREE, bond_distribution_ALLFREELIMSG, ]

        if target_data[index]["bond_detail_level"] == "Complete":
            axes = ax[0]
            ax0labels.append(biomass_labels[index])
            ax0index += 1
            barindex = ax0index


        if target_data[index]["bond_detail_level"] == "terrell":
            axes = ax[1]
            ax1labels.append(biomass_labels[index])
            ax1index += 1
            barindex = ax1index

        if target_data[index]["bond_detail_level"] == "klose":
            axes = ax[2]
            ax2labels.append(biomass_labels[index])
            ax2index +=1
            barindex = ax2index
        bottomtarget = 0
        bottom = np.zeros(3)
        fitness = np.zeros(3)

        for j in range(len(bond_dist_list)):
            bond_distribution = bond_dist_list[j]
            BO4[j] = bond_distribution["bo4"]
            BB[j] = bond_distribution["bb"]
            B5[j] = bond_distribution["b5"]
            B1[j] = bond_distribution["b1"]
            AO4[j] = bond_distribution["ao4"]
            _4O5[j] = bond_distribution["5o4"]
            _55[j] = bond_distribution["55"]
            OTHER[j] = bond_distribution["ao4"] + bond_distribution["5o4"] + bond_distribution["55"]
            fitness[j] = round(bond_metrics(target_data[index], bond_distribution)["euclidian"], 2)

            if target_data[index]["bond_detail_level"] == "klose":
                bondsum = bond_distribution["bo4"] +  bond_distribution["bb"] + bond_distribution["b5"]
                BO4[j] /= bondsum
                BB[j] /= bondsum
                B5[j] /= bondsum
        axes.bar(barindex-1, target_BO4, barwidth - subtr , edgecolor='white', color=colors[0])
        bottomtarget += target_BO4
        axes.bar(barindex-1 , target_BB, barwidth - subtr, bottom=bottomtarget,  edgecolor='white', color=colors[1])
        bottomtarget += target_BB
        axes.bar(barindex-1 , target_B5, barwidth - subtr, bottom=bottomtarget,edgecolor='white', color=colors[2])
        bottomtarget += target_B5
        axes.text(barindex-1 , 1.02, "Plant", fontsize=12, horizontalalignment = 'center')
        #print(barindex, BO4[index])
        for sim in range(len(library_list)):
            axes.text(barindex-1+ separation*sim + barwidth + dist, 1.02, str(fitness[sim]), fontsize=12,horizontalalignment = 'center')
            axes.text(barindex-1+ separation*sim + barwidth + dist, -0.1, plotlablist[sim], fontsize=15,horizontalalignment = 'center')
            axes.bar(barindex-1+ separation*sim + barwidth + dist, BO4[sim], barwidth,  edgecolor='white', color=colors[0])
            bottom[sim] = BO4[sim]
            axes.bar(barindex-1+ separation*sim + barwidth + dist, BB[sim], barwidth, bottom=bottom[sim], edgecolor='white', color=colors[1])
            bottom[sim] += BB[sim]
            axes.bar(barindex-1+ separation*sim + barwidth + dist, B5[sim], barwidth, bottom=bottom[sim], edgecolor='white', color=colors[2])
            bottom[sim] += B5[sim]
        
        if target_data[index]["bond_detail_level"] == "Complete":
            for sim in range(len(library_list)):
                axes.bar(barindex-1+ separation*sim + barwidth + dist, _55[sim], barwidth, bottom=bottom[sim], edgecolor='white', color=colors[3])
                bottom[sim] += _55[sim]
                axes.bar(barindex-1+ separation*sim + barwidth + dist, AO4[sim], barwidth, bottom=bottom[sim],  edgecolor='white', color=colors[4])
                bottom[sim] += AO4[sim]
                axes.bar(barindex-1+ separation*sim + barwidth + dist, _4O5[sim], barwidth, bottom=bottom[sim], edgecolor='white', color=colors[5])
                bottom[sim] += _4O5[sim]


            axes.bar(barindex-1 , target_55, barwidth - subtr, bottom=bottomtarget, edgecolor='white', color=colors[3])
            bottomtarget += target_55
            axes.bar(barindex-1 , target_AO4, barwidth - subtr, bottom=bottomtarget, edgecolor='white', color=colors[4])
            bottomtarget += target_AO4
            axes.bar(barindex-1 , target_4O5, barwidth - subtr, bottom=bottomtarget,  edgecolor='white', color=colors[5])
            bottomtarget += target_4O5
                

    
        if target_data[index]["bond_detail_level"] == "Complete" or target_data[index]["bond_detail_level"] == "terrell":
            for sim in range(len(library_list)):
                axes.bar(barindex-1+ separation*sim + barwidth + dist, B1[sim], barwidth, bottom=bottom[sim], edgecolor='white', color=colors[6])
                bottom[sim] += B1[sim]

            axes.bar(barindex-1 , target_B1, barwidth - subtr, bottom=bottomtarget,  edgecolor='white', color=colors[6])
            bottomtarget += target_B1

        if target_data[index]["bond_detail_level"] == "terrell":
            for sim in range(len(library_list)):
                axes.bar(barindex-1+ separation*sim + barwidth + dist, OTHER[sim], barwidth, bottom=bottom[sim], edgecolor='white', color='lightgrey')
                bottom[sim] += OTHER[sim]
            axes.bar(barindex-1 , target_OTHER, barwidth - subtr, bottom=bottomtarget, edgecolor='white', color='lightgrey')
            bottomtarget += target_OTHER

        index += 1
    ##PLOT LEGEND
    axes = ax[1]
    axes.axvline(barindex-1-0.15, color="k", linestyle="--", linewidth=2)
    ax1labels.append("ID")
    ax1index +=1
    barindex = ax1index
    sample = np.array([0.4, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
    bottom = np.zeros(3)
    bottomtarget = 0
    axes.text(barindex-1+0.25, 1.02, "Experiment", fontsize=12, horizontalalignment = 'center')
    axes.text(barindex-1 +0.7, 1.02, "Fitness", fontsize=12, horizontalalignment = 'center')
    axes.text(barindex-1 +0.7, -0.15, "Simulated\nDistributions", fontsize=12, horizontalalignment = 'center')

    for num, test in enumerate(sample):
        axes.bar(barindex-1 +0.25 , test, barwidth - subtr, bottom=bottomtarget, edgecolor='black', color=colors[num])
        bottomtarget += test
        for i in range(len(library_list)):
            axes.bar(barindex-1+ 0.25+separation*i + barwidth + dist, test, barwidth, bottom=bottom[sim], alpha=0.7, edgecolor='black', color=colors[num])
            bottom[i] += sample[num]



    pos =0#(width - dist)
    ax[0].set_xticks(np.arange(pos, ax0index+pos))
    ax[0].set_xticklabels(ax0labels, fontsize=fontsize+2, fontdict={'weight': 'bold'})#, rotation=30)
    ax[2].set_xticks(np.arange(pos, ax2index+pos))
    ax[2].set_xticklabels(ax2labels, fontsize=fontsize+2, fontdict={'weight': 'bold'})# rotation=30)
    arr = np.arange(pos, ax1index-1+pos)
    arr = np.append(arr, ax1index-1+0.25)
    ax[1].set_xticks(arr)
    ax[1].set_xticklabels(ax1labels, fontsize=fontsize+2, fontdict={'weight': 'bold'})# rotation=30)

    ax[0].set_yticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], fontsize=fontsize-2)# rotation=30)
    ax[1].set_yticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], fontsize=fontsize-2)# rotation=30)
    ax[2].set_yticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], fontsize=fontsize-2)# rotation=30)
    ax[0].set_ylim(0, 1.4)
    ax[1].set_ylim(0, 1.4)
    ax[2].set_ylim(0, 1.4)
    # Define the legend labels and colors
    legend_labels = ["$\\beta-O-4$", "$\\beta-\\beta$", "$\\beta-5$", "$5-5$", "$\\alpha-O-4$", "$4-O-5$", "$\\beta-1$"]
    legend_colors = colors[:len(legend_labels)]

    # Create custom legend handles using Patch for stacked bars
    legend_handles = [patches.Patch(color=color, label=label) for color, label in zip(legend_colors, legend_labels)]
    legend_handles.append(patches.Patch(color="lightgrey", label="Sum of $\\alpha-O-4$, $5-5$\nand $4-O-5$"))
    legend_labels.append("Sum of $\\alpha-O-4$, $5-5$\nand $4-O-5$")
    legend_handles.append(patches.Patch(color="white", label="Simulation Conditions"))
    legend_labels.append("Simulation Conditions")    
    legend_handles.append(patches.Patch(color="white", label="Simulation Conditions"))
    legend_labels.append("A: Monomer addition rate and S/G ratio free")    
    legend_handles.append(patches.Patch(color="white", label="A: Monomer addition rate and S/G ratio free"))
    legend_labels.append("B: All parameters free")    
    legend_handles.append(patches.Patch(color="white", label="B: All parameters free"))
    legend_labels.append("C: All parameters free, limited S/G")
    axes.plot([],[], label ="")
    axes.plot([],[], label ="")
    axes.plot([],[], label="")


    #bond_distribution_MONADDSG, bond_distribution_ALLFREE, bond_distribution_ALLFREELIMSG,
    # Create a legend for the main plot
    legend = fig.legend(handles=legend_handles, labels=legend_labels, loc='lower center', bbox_to_anchor=(0.5, -0.2), fontsize=fontsize, ncol=3)

    # Add the second legend to the plot
    #extra = plt.gca().add_artist(legend1)
    ax[0].set_ylim(0, 1.1)
    savestring = "bond_distributions_and_targets.png"
        # Adjust the layout to prevent overlapping titles
    plt.tight_layout()
    plt.savefig(savestring, bbox_inches='tight')





def plot_compare_euclidian_distances(libraries_list, biomasses, labels):
    """
    Compares the old and new euclidian distances (y axis) of different fits for N different biomasses (x axis)
    """
    old_libraries = libraries_list[0]
    mid_libraries = libraries_list[1]
    new_libraries = libraries_list[2]

    fontsize = 22
    fig, ax = plt.subplots(figsize=(12,8))
     # Set the axis labels and title
    ax.set_xlabel('Biomass', fontsize=fontsize)
    ax.set_ylabel('Euclidian distance from target', fontsize=fontsize)

    colors = [cmap(i) for i in np.linspace(0, 1, 10)]
    index=0
    color1 = colors[0]
    color2 = colors[3]
    color3 = colors[5]
    color4 = colors[7]
    # Set positions of each biomass
    x_positions = np.arange(len(new_libraries))
    ax.set_xticks(x_positions)
    ax.set_xticklabels(labels)
    ax.tick_params(axis='x', labelsize=fontsize)  # Change the fontsize of x-axis ticks
    ax.tick_params(axis='y', labelsize=fontsize)  # Change the fontsize of y-axis ticks
    #new_euclidian_distances =  np.zeros(len(new_libraries))
    #old_euclidian_distances =  np.zeros(len(old_libraries))
    transparency = 0.8

    for ind, library in new_libraries:
        # calculate average bond distributon
        simulated_distribution = calculate_avg_bond_distribution(library)
       
        # calculate euclidian distance to biomass
        new_euclidian_distance =  bond_metrics(biomasses[index], simulated_distribution)["euclidian"]
        if index == 0:
            ax.scatter(x_positions[index], new_euclidian_distance, marker="^", alpha=transparency, color = color1, s=300, label = "Monomer addition rate, SG ratio free")

        else:
            ax.scatter(x_positions[index], new_euclidian_distance, marker="^", alpha=transparency, color = color1, s=300)
        index+=1

    index=0
    for ind, library in old_libraries:
        # calculate average bond distributon
        simulated_distribution = calculate_avg_bond_distribution(library)
        # calculate euclidian distance to biomass
        old_euclidian_distance =  bond_metrics(biomasses[index], simulated_distribution)["euclidian"]
        if index == 0:
            ax.scatter(x_positions[index], old_euclidian_distance, marker="^", alpha=transparency, color = color2, s=300, label = "All parameters free")

        else:
            ax.scatter(x_positions[index], old_euclidian_distance, marker="^", alpha=transparency, color = color2, s=300)
        index+=1
    index=0
    for ind, library in mid_libraries:
        # calculate average bond distributon
        simulated_distribution = calculate_avg_bond_distribution(library)
        # calculate euclidian distance to biomass
        old_euclidian_distance =  bond_metrics(biomasses[index], simulated_distribution)["euclidian"]
        if index == 0:
            ax.scatter(x_positions[index], old_euclidian_distance, marker="^", alpha=transparency, color = color3, s=300, label = "All parameters free, limited SG")

        else:
            ax.scatter(x_positions[index], old_euclidian_distance, marker="^", alpha=transparency, color = color3, s=300)
        index+=1

    savestring = "compare_fit_quality.png"
    # Adjust the layout to prevent overlapping titles    
    fig.legend(fontsize=18, loc='lower center', bbox_to_anchor=(0.5, -0.12), markerscale = 1.5, ncols=2)

    plt.tight_layout()
    plt.savefig(savestring, bbox_inches='tight')

def calculate_weight_avg(library):
    """
    Calculates the weight-average molecular weight (Mw) of a library of molecules.

    Parameters:
    -----------
    library : list of dicts
        A list of dictionaries where each dictionary contains information about a molecule,
        including a 'smilestring' key that provides the SMILES representation of the molecule.

    Returns:
    --------
    MW : float
        The weight-average molecular weight (Mw) of the molecules in the library.

    Description:
    ------------
    This function computes the weight-average molecular weight (Mw) for a given library of molecules.
    It extracts the SMILES strings from the library and calculates the molecular weight for each molecule.
    The molecular weights are then counted to determine the frequency of each unique molecular weight.

    The function calculates the weight fractions for each molecular weight, which are used to compute
    two quantities: the total weight (sum of molecular weights multiplied by their fractions) and
    the total weight times the molecular weight (sum of the square of molecular weights multiplied by their fractions).

    The weight-average molecular weight (Mw) is then calculated as the ratio of the total weight times the molecular weight
    to the total weight. This value is returned as the output of the function.
    """

    smiles = [data["smilestring"] for data in library]
    molecular_weights = [Descriptors.MolWt(MolFromSmiles(smilestring)) for smilestring in smiles]


    weight_counts = Counter(molecular_weights)
    
    # Total number of molecules
    total_count = sum(weight_counts.values())
    
    # Calculate the fractions (weights)
    weight_fractions = {mw: count / total_count for mw, count in weight_counts.items()}
    
    # Calculate total_weight and total_weight_times_fraction
    total_weight = sum(mw * fraction for mw, fraction in weight_fractions.items())
    total_weight_times_fraction = sum(mw**2 * fraction for mw, fraction in weight_fractions.items())

    # Weight average molecular weight (Mw)
    MW = total_weight_times_fraction / total_weight
    
    return(MW)

def calculate_number_avg(library):
    """
    Calculates the number-average molecular weight (Mn) of a library of molecules.

    Parameters:
    -----------
    library : list of dicts
        A list of dictionaries where each dictionary contains information about a molecule,
        including a 'smilestring' key that provides the SMILES representation of the molecule.

    Returns:
    --------
    Mn : float
        The number-average molecular weight (Mn) of the molecules in the library.

    Description:
    ------------
    This function computes the number-average molecular weight (Mn) for a given library of molecules.
    It first extracts the SMILES strings from the library and calculates the molecular weight for each molecule.
    The molecular weights are then counted to determine the frequency of each unique molecular weight.
    
    The function computes the total number of molecules and uses these counts to calculate the weight fractions
    for each molecular weight. Finally, the function returns Mn, which is calculated as the sum of the molecular
    weights multiplied by their respective fractions.
    """
    smiles = [data["smilestring"] for data in library]
    molecular_weights = [Descriptors.MolWt(MolFromSmiles(smilestring)) for smilestring in smiles]
    weight_counts = Counter(molecular_weights)
    
    # Total number of molecules
    total_count = sum(weight_counts.values())
    
    # Calculate the fractions (weights)
    weight_fractions = {mw: count / total_count for mw, count in weight_counts.items()}
    
    # Calculate total_weight and total_weight_times_fraction
    total_weight = sum(mw * fraction for mw, fraction in weight_fractions.items())

    # Weight average molecular weight (Mw)
    Mn= total_weight
    return(Mn)

def plot_polydispersityindex_comparison(library_list, biomasses, labels):
    """
    Plots a comparison of molecular weight averages ($M_w$), number averages ($M_n$), 
    and polydispersity indices (PDI) across different simulation conditions and experimental data.

    Parameters:
    -----------
    library_list : list of lists
        A list containing sublists, where each sublist corresponds to a different simulation condition.
        Each sublist contains data needed to calculate $M_w$, $M_n$, and PDI for different biomasses.
    biomasses : list of dicts
        A list of dictionaries where each dictionary contains experimental data for a particular biomass.
        The dictionary may contain keys 'MW' for molecular weight, 'MN' for number average, and 'Monolignols' for specific monomer ratios.
    labels : list of str
        A list of labels for each biomass, used for the x-axis of the plot.

    Returns:
    --------
    None
        The function generates and saves a plot 'compare_PDI.png' that displays $M_w$, $M_n$, and PDI 
        for each biomass and simulation condition. The plot includes experimental data for comparison.

    Description:
    ------------
    This function creates a figure with three subplots, each corresponding to a different metric:
    1. $M_w$ (Weight-average molecular weight)
    2. $M_n$ (Number-average molecular weight)
    3. PDI (Polydispersity Index), calculated as $M_w / M_n$

    For each biomass in `biomasses`, the function computes these metrics across multiple simulation conditions 
    provided in `library_list`. Experimental data (if available) is also plotted for comparison.

    The subplots are aligned vertically, share the same x-axis, and are labeled with the corresponding biomass labels. 
    The plot includes a custom legend that describes the different simulation conditions and experimental data.

    The plot is customized with colors, labels, and legend, and is saved as 'compare_PDI.png'.
    """
    fontsize = 22
    fig, axes = plt.subplots(3,figsize=(15,8), sharex=True)
     # Set the axis labels and title
    axes[2].set_xlabel('Biomass', fontsize=fontsize)
    axes[0].set_ylabel('$M_w$', fontsize=fontsize)
    axes[1].set_ylabel('$M_n$', fontsize=fontsize)
    axes[2].set_ylabel('PDI', fontsize=fontsize)

    # Bar position and shape properties 
    dist = 0.0
    subtr = -0.05
    separation = 0.2
    width = 0.65  # Reduce the width for the target distribution
    barwidth = 0.22
    cmap = plt.get_cmap('gnuplot2')
    colors = [cmap(i) for i in np.linspace(0, 1, 12)]
    color1 = colors[4]
    color2 = colors[6]
    color3 = colors[8]
    barindex=0
    index = 0
    S_monolignol_weight = 184.19
    G_monolignol_weight = 124.14
    color4 = colors[2]
    colorlist= [color1, color2, color3]
    # Set positions of each biomass
    x_positions = np.arange(len(library_list[0]))
    axes[2].set_xticks(x_positions-1 + 0.25)
    axes[2].set_xticklabels(labels)

    plotlablist = ["A", "B", "C"]

    #new_euclidian_distances =  np.zeros(len(new_libraries))
    #old_euclidian_distances =  np.zeros(len(old_libraries))
    transparency = 1
    for i in range(len(library_list[0])):
        print(i)
        print("All params free")
        mn_ALLFREE = calculate_number_avg(library_list[0][index][1][:500])
        mw_ALLFREE = calculate_weight_avg(library_list[0][index][1][:500])
        pd_ALLFREE = mw_ALLFREE/mn_ALLFREE
        #pd_ALLFREE = np.average(DPs)#, weights = DPs)
        print("All params free limsg")

        mn_ALLFREELIMSG = calculate_number_avg( library_list[1][index][1][:500])
        mw_ALLFREELIMSG = calculate_weight_avg(library_list[1][index][1][:500])
        pd_ALLFREELIMSG = mw_ALLFREELIMSG/mn_ALLFREELIMSG

        print("MonaddSG")
        mn_MONADDSG = calculate_number_avg(library_list[2][index][1][:500])
        mw_MONADDSG = calculate_weight_avg(library_list[2][index][1][:500])
        pd_MONADDSG = mw_MONADDSG/mn_MONADDSG


        mw_list = [mw_MONADDSG, mw_ALLFREE, mw_ALLFREELIMSG]
        mn_list = [mn_MONADDSG, mn_ALLFREE, mn_ALLFREELIMSG]
        pd_list = [pd_MONADDSG, pd_ALLFREE, pd_ALLFREELIMSG]
        for sim in range(len(library_list)):
            #fitness[sim] = round(bond_metrics(target_data[index], bond_distribution)["euclidian"], 2)
            #axes.text(barindex-1+ separation*sim + barwidth + dist, sg_list[sim] + 0.1, str(fitness[sim]), rotation=90, fontsize=11,horizontalalignment = 'center')
            #axes[0].text(barindex-1+ separation*sim + barwidth + dist, -0.25, plotlablist[sim], fontsize=12,horizontalalignment = 'center')
            axes[0].bar(barindex-1+ separation*sim + barwidth + dist, mw_list[sim], barwidth,  edgecolor='white', color=colorlist[sim])
            #axes[1].text(barindex-1+ separation*sim + barwidth + dist, -0.25, plotlablist[sim], fontsize=12,horizontalalignment = 'center')
            axes[1].bar(barindex-1+ separation*sim + barwidth + dist, mn_list[sim], barwidth,  edgecolor='white', color=colorlist[sim])
            #axes[2].text(barindex-1+ separation*sim + barwidth + dist, -0.25, plotlablist[sim], fontsize=12,horizontalalignment = 'center')
            axes[2].bar(barindex-1+ separation*sim + barwidth + dist, pd_list[sim], barwidth,  edgecolor='white', color=colorlist[sim])
        if 'MN' in biomasses[index]:
            #target_DP = (biomasses[index]["MN"]*(biomasses[index]["Monolignols"]["S"]/100)/S_monolignol_weight) + (biomasses[index]["MN"]*(biomasses[index]["Monolignols"]["G"]/100)/G_monolignol_weight)
            target = biomasses[index]["MN"]
            axes[1].bar(barindex-1 , target, barwidth - subtr,  edgecolor='white', color="k")
        if 'MW' in biomasses[index]:
            #target_DP = (biomasses[index]["MN"]*(biomasses[index]["Monolignols"]["S"]/100)/S_monolignol_weight) + (biomasses[index]["MN"]*(biomasses[index]["Monolignols"]["G"]/100)/G_monolignol_weight)
            target = biomasses[index]["MW"]
            axes[0].bar(barindex-1 , target, barwidth - subtr,  edgecolor='white', color="k")
        if 'MN' in biomasses[index] and 'MW' in biomasses[index]:
            #target_DP = (biomasses[index]["MN"]*(biomasses[index]["Monolignols"]["S"]/100)/S_monolignol_weight) + (biomasses[index]["MN"]*(biomasses[index]["Monolignols"]["G"]/100)/G_monolignol_weight)
            target = biomasses[index]["MW"]/biomasses[index]["MN"]
            axes[2].bar(barindex-1 , target, barwidth - subtr,  edgecolor='white', color="k")
        index +=1
        barindex +=1


    # Define the legend labels and colors
    legend_labels = ["Experimental data","A: Monomer addition rate and S/G ratio free", "B: All parameters free", "C: All parameters free, limited S/G"]
    legend_colors = ["k", color1, color2, color3]

    # Create custom legend handles using Patch for stacked bars
    legend_handles = [patches.Patch(color=color, label=label) for color, label in zip(legend_colors, legend_labels)]



    #bond_distribution_MONADDSG, bond_distribution_ALLFREE, bond_distribution_ALLFREELIMSG,
    # Create a legend for the main plot
    legend = fig.legend(handles=legend_handles, labels=legend_labels, loc='lower center', bbox_to_anchor=(0.55, -0.2), fontsize=fontsize, ncol=2)
    # Adjust the layout to prevent overlapping titles    
    axes[2].tick_params(axis='x', labelsize=fontsize)  # Change the fontsize of x-axis ticks
    axes[0].tick_params(axis='y', labelsize=fontsize)  # Change the fontsize of y-axis ticks
    axes[1].tick_params(axis='y', labelsize=fontsize)  # Change the fontsize of y-axis ticks    axes[0].tick_params(axis='x', labelsize=fontsize)  # Change the fontsize of x-axis ticks
    axes[2].tick_params(axis='y', labelsize=fontsize)  # Change the fontsize of y-axis ticks
        

    plt.tight_layout()
    savestring = "compare_PDI.png"
    plt.savefig(savestring, bbox_inches='tight')
    return()

#def zimmstockmeyer_tribranching_heterogeneous(n_branches_weighted):
 #   nw = n_branches_weighted
  #  term_1 = 0.5*(np.sqrt((2+nw)/(nw)))
   # term_2_numerator = (np.sqrt(2+nw)) + (np.sqrt(nw))
   # term_2_denominator = (np.sqrt(2+nw)) - (np.sqrt(nw))
  #  term_2 = (np.log((term_2_numerator/term_2_denominator))) - 1
   # g3w = (term_1*term_2) * (6/nw)

    #return(g3w)

def zimmstockmeyer_tribranching_homogeneous(branches):
    b = branches
    g3 = ( np.sqrt(1+(b/7)) + (4*b)/(9*np.pi) )**(-0.5)
    return(g3)

def plot_branching_subplots(library_list, biomasses, labels):

    fontsize = 22
    fig, axes = plt.subplots(3,figsize=(15,8), sharex=True)
     # Set the axis labels and title
    axes[2].set_xlabel('Biomass', fontsize=fontsize)
    axes[0].set_ylabel('Branching\ncoefficient $\\alpha$', fontsize=fontsize)
    axes[1].set_ylabel('DP weighted\nbranching factor', fontsize=fontsize)
    axes[2].set_ylabel('Branching\nratio $g_{3}$', fontsize=fontsize)

    # Bar position and shape properties 
    dist = 0.0
    subtr = -0.05
    separation = 0.2
    width = 0.65  # Reduce the width for the target distribution
    barwidth = 0.22
    cmap = plt.get_cmap('gnuplot2')
    colors = [cmap(i) for i in np.linspace(0, 1, 12)]
    color1 = colors[4]
    color2 = colors[6]
    color3 = colors[8]
    barindex=1
    index = 0

    color4 = colors[2]
    colorlist= ["k", color1, color2, color3]
    # Set positions of each biomass
    x_positions = np.arange(1, len(library_list[0])+1)
    axes[2].set_xticks(x_positions + 0.25)
    axes[2].set_xticklabels(labels)

    plotlablist = ["A", "B", "C"]

    #new_euclidian_distances =  np.zeros(len(new_libraries))
    #old_euclidian_distances =  np.zeros(len(old_libraries))
    transparency = 1
    branch = []
    zimmstockmeyer_branching = []
    alphas =[]
    DPs = []
    cmap = plt.get_cmap('rainbow')

    colors = [cmap(i) for i in np.linspace(0, 1, len(library_list[0]))]
    for i in range(len(library_list[0])):
        print(i)
        mn_ALLFREELIMSG = calculate_number_avg(library_list[1][index][1][:5])
        mw_ALLFREELIMSG = calculate_weight_avg(library_list[1][index][1][:5])
        pd_ALLFREELIMSG = mw_ALLFREELIMSG/mn_ALLFREELIMSG
        pd_inv = 1/pd_ALLFREELIMSG#[pd_MONADDSG, pd_ALLFREE, pd_ALLFREELIMSG]
        if pd_inv < 0.63:
            alpha = (-0.58*pd_inv) + 0.36
        else:
            alpha = 0 
        alphas.append(alpha)
        #alphas.append(0)
        #alpha =0


            #axes.text(barindex-1+ separation*sim + barwidth + dist, sg_list[sim] + 0.1, str(fitness[sim]), rotation=90, fontsize=11,horizontalalignment = 'center')
            #axes[0].text(barindex-1+ separation*sim + barwidth + dist, -0.25, plotlablist[sim], fontsize=12,horizontalalignment = 'center')
        axes[0].bar(barindex + barwidth + dist, alpha, barwidth,  edgecolor='white', color=colorlist[3])
        #if 'MN' in biomasses[index] and 'MW' in biomasses[index]:

        if 'alpha' in biomasses[index]:
            #target_DP = (biomasses[index]["MN"]*(biomasses[index]["Monolignols"]["S"]/100)/S_monolignol_weight) + (biomasses[index]["MN"]*(biomasses[index]["Monolignols"]["G"]/100)/G_monolignol_weight)
            target = biomasses[index]["alpha"]#["MN"]/biomasses[index]["MW"]
            if target <0.05:
                target = 0.01
            axes[0].bar(barindex , target, barwidth - subtr,  edgecolor='white', color="k")

        elif "MW" in biomasses[index] and "MN" in biomasses[index]:
            pd_inv = biomasses[index]["MN"]/biomasses[index]["MW"]#[pd_MONADDSG, pd_ALLFREE, pd_ALLFREELIMSG]
            if pd_inv < 0.63:
                target = (-0.58*pd_inv) + 0.36
            else:
                alpha = 0 
            if target <0.05:
                target = 0.01
            axes[0].bar(barindex , target, barwidth - subtr,  edgecolor='white', color="k")
        _branching = np.array([data["Branches"] for data in library_list[1][index][1]])
        DP = np.array([data["DP"] for data in library_list[1][index][1]])
        # Limit branching to the first 10 elements
        #branching = _branching[:400]

        # Extract the first 10 SMILES strings
        """smiles = [data["smilestring"] for data in library_list[1][index][1][:400]]

        # Calculate molecular weights
        molecular_weights = [Descriptors.MolWt(MolFromSmiles(smilestring)) for smilestring in smiles]

        # Count the occurrences of each molecular weight and maintain the order
        weight_counts = OrderedDict()
        for mw in molecular_weights:
            if mw in weight_counts:
                weight_counts[mw] += 1
            else:
                weight_counts[mw] = 1

        # Total number of molecules
        total_count = sum(weight_counts.values())

        # Calculate the fractions (weights) of each molecular weight
        #weight_fractions = OrderedDict((mw, count / total_count) for mw, count in weight_counts.items())

        # Calculate the weighted branching values and sum them up
        weighted_sum = 0
        #branches = np.zeros(len(weight_fractions.items()))
        #for ind, (mw, fraction) in enumerate(weight_fractions.items()):
        #    weighted_sum += branching[ind] * fraction
        #    branches[ind] = branching[ind] * fraction

        # The weight average number of branches
        #weight_average_branches = weighted_sum
        #print(weight_average_branches)
        #print(branches)"""
        zimmstockmeyer_branch = zimmstockmeyer_tribranching_homogeneous(_branching)#_weight_average_branches)
        zimmstockmeyer_branching.append(zimmstockmeyer_branch)
        DPs.append(DP)
        # Create weighted data by repeating each data point according to its weight
        weighted_data = np.repeat(_branching/DP, DP)
        branch.append(weighted_data)
        #axes[2].bar(barindex + barwidth + dist, zimmstockmeyer_branch, barwidth,  edgecolor='white', color=colors[i])
        #axes[2].set_ylim(0, 2)

        index+=1

        barindex +=1

    violin_parts = axes[1].violinplot(branch,showmeans=True, showextrema=True)#
    violin_parts_2 = axes[2].violinplot(zimmstockmeyer_branching,showmeans=True, showextrema=True)#


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
    if 'cmeans' in violin_parts:
        violin_parts['cmeans'].set_color('black')
    # Customize colors

    for i, pc in enumerate(violin_parts_2['bodies']):
        pc.set_facecolor("grey")
        pc.set_edgecolor("grey")
        pc.set_alpha(0.7)


    # Customize range lines (extrema) if they exist
    if 'cmins' in violin_parts_2:
        violin_parts_2['cmins'].set_color('black')
    if 'cmaxes' in violin_parts:
        violin_parts_2['cmaxes'].set_color('black')
    if 'cbars' in violin_parts:
        violin_parts_2['cbars'].set_color('black')
    if 'cmeans' in violin_parts:
        violin_parts_2['cmeans'].set_color('black')


    ####### FOR THE ALPHAS
    means_alpha = np.array(alphas)
    data_alpha = {
        'Group': ['Grasses'] * 4 + ['Hardwoods'] * 7 + ['Softwoods'] * 5,
        'BranchingParameter': means_alpha[:16]
    }

    df_alpha = pd.DataFrame(data_alpha)
# Perform ANOVA
    anova_result = stats.f_oneway(
        df_alpha[df_alpha['Group'] == 'Grasses']['BranchingParameter'],
        df_alpha[df_alpha['Group'] == 'Hardwoods']['BranchingParameter'],
        df_alpha[df_alpha['Group'] == 'Softwoods']['BranchingParameter']
    ) 
    print(f"F-statistic: {anova_result.statistic}, p-value: {anova_result.pvalue}")
    # Perform Tukey's HSD test
    tukey_result = pairwise_tukeyhsd(df_alpha['BranchingParameter'], df_alpha['Group'], alpha=0.05)
    print("Test results for alpha")
    print(tukey_result) 


    ####### FOR THE BRANCHING FACTOR
    means_bf = [segment[0, 1] for segment in violin_parts['cmeans'].get_segments()]

    data_bf = {
        'Group': ['Grasses'] * 4 + ['Hardwoods'] * 7 + ['Softwoods'] * 5,
        'BranchingParameter': means_bf[:16]
    }

    df_bf = pd.DataFrame(data_bf)
# Perform ANOVA
    anova_result = stats.f_oneway(
        df_bf[df_bf['Group'] == 'Grasses']['BranchingParameter'],
        df_bf[df_bf['Group'] == 'Hardwoods']['BranchingParameter'],
        df_bf[df_bf['Group'] == 'Softwoods']['BranchingParameter']
    ) 
    print(f"F-statistic: {anova_result.statistic}, p-value: {anova_result.pvalue}")
    # Perform Tukey's HSD test
    tukey_result = pairwise_tukeyhsd(df_bf['BranchingParameter'], df_bf['Group'], alpha=0.05)
    print("Test results for branching factor")
    print(tukey_result) 

    ####### FOR THE G3
    means_g3 = [segment[0, 1] for segment in violin_parts_2['cmeans'].get_segments()]
    data_g3 = {
        'Group': ['Grasses'] * 4 + ['Hardwoods'] * 7 + ['Softwoods'] * 5,
        'BranchingParameter': means_g3[:16]
    }

    df_g3 = pd.DataFrame(data_g3)
# Perform ANOVA
    anova_result = stats.f_oneway(
        df_g3[df_g3['Group'] == 'Grasses']['BranchingParameter'],
        df_g3[df_g3['Group'] == 'Hardwoods']['BranchingParameter'],
        df_g3[df_g3['Group'] == 'Softwoods']['BranchingParameter']
    ) 
    print(f"F-statistic: {anova_result.statistic}, p-value: {anova_result.pvalue}")
    # Perform Tukey's HSD test
    tukey_result = pairwise_tukeyhsd(df_g3['BranchingParameter'], df_g3['Group'], alpha=0.05)
    print("Test results for g3")
    print(tukey_result) 
    

    # Define the legend labels and colors
    legend_labels = ["Experimental data", "C: All parameters free, limited S/G"]
    legend_colors = ["k", color3]

    # Create custom legend handles using Patch for stacked bars
    legend_handles = [patches.Patch(color=color, label=label) for color, label in zip(legend_colors, legend_labels)]

    #bond_distribution_MONADDSG, bond_distribution_ALLFREE, bond_distribution_ALLFREELIMSG,
    # Create a legend for the main plot
    legend = fig.legend(handles=legend_handles, labels=legend_labels, loc='lower center', bbox_to_anchor=(0.6, -0.2), fontsize=fontsize, ncol=2)
    # Adjust the layout to prevent overlapping titles    
    axes[2].tick_params(axis='x', labelsize=fontsize)  # Change the fontsize of x-axis ticks
    axes[0].tick_params(axis='y', labelsize=fontsize)  # Change the fontsize of y-axis ticks
    axes[1].tick_params(axis='y', labelsize=fontsize)  # Change the fontsize of y-axis ticks    axes[0].tick_params(axis='x', labelsize=fontsize)  # Change the fontsize of x-axis ticks
    axes[2].tick_params(axis='y', labelsize=fontsize)  # Change the fontsize of y-axis ticks
        

    plt.tight_layout()
    savestring = "compare_branching_tri.png"
    plt.savefig(savestring, bbox_inches='tight')
    return()


def plot_g3_branchingfactor_subplots(library_list, biomasses, labels):

    fontsize = 22
    fig, axes = plt.subplots(2,figsize=(14,7), sharex=True)
     # Set the axis labels and title
    axes[1].set_xlabel('Biomass', fontsize=fontsize)
    #axes[0].set_ylabel('Branching\ncoefficient $\\alpha$', fontsize=fontsize)
    axes[0].set_ylabel('DP weighted\nbranching factor', fontsize=fontsize)
    axes[1].set_ylabel('Branching\nratio $g_{3}$', fontsize=fontsize)

    # Bar position and shape properties 
    dist = 0.0
    subtr = -0.05
    separation = 0.2
    width = 0.65  # Reduce the width for the target distribution
    barwidth = 0.22
    cmap = plt.get_cmap('gnuplot2')
    colors = [cmap(i) for i in np.linspace(0, 1, 12)]
    color1 = colors[4]
    color2 = colors[6]
    color3 = colors[8]
    barindex=1
    index = 0

    color4 = colors[2]
    colorlist= ["k", color1, color2, color3]
    # Set positions of each biomass
    x_positions = np.arange(1, len(library_list[0])+1)
    axes[1].set_xticks(x_positions)
    axes[1].set_xticklabels(labels)

    plotlablist = ["A", "B", "C"]

    #new_euclidian_distances =  np.zeros(len(new_libraries))
    #old_euclidian_distances =  np.zeros(len(old_libraries))
    transparency = 1
    branch = []
    zimmstockmeyer_branching = []
    alphas =[]
    DPs = []
    cmap = plt.get_cmap('rainbow')

    colors = [cmap(i) for i in np.linspace(0, 1, len(library_list[0]))]
    for i in range(len(library_list[0])):
        print(i)


        _branching = np.array([data["Branches"] for data in library_list[1][index][1]])
        DP = np.array([data["DP"] for data in library_list[1][index][1]])
        # Limit branching to the first 10 elements
        #branching = _branching[:400]

        # Extract the first 10 SMILES strings
        """smiles = [data["smilestring"] for data in library_list[1][index][1][:400]]

        # Calculate molecular weights
        molecular_weights = [Descriptors.MolWt(MolFromSmiles(smilestring)) for smilestring in smiles]

        # Count the occurrences of each molecular weight and maintain the order
        weight_counts = OrderedDict()
        for mw in molecular_weights:
            if mw in weight_counts:
                weight_counts[mw] += 1
            else:
                weight_counts[mw] = 1

        # Total number of molecules
        total_count = sum(weight_counts.values())

        # Calculate the fractions (weights) of each molecular weight
        #weight_fractions = OrderedDict((mw, count / total_count) for mw, count in weight_counts.items())

        # Calculate the weighted branching values and sum them up
        weighted_sum = 0
        #branches = np.zeros(len(weight_fractions.items()))
        #for ind, (mw, fraction) in enumerate(weight_fractions.items()):
        #    weighted_sum += branching[ind] * fraction
        #    branches[ind] = branching[ind] * fraction

        # The weight average number of branches
        #weight_average_branches = weighted_sum
        #print(weight_average_branches)
        #print(branches)"""
        zimmstockmeyer_branch = zimmstockmeyer_tribranching_homogeneous(_branching)#_weight_average_branches)
        zimmstockmeyer_branching.append(zimmstockmeyer_branch)
        DPs.append(DP)
        # Create weighted data by repeating each data point according to its weight
        weighted_data = np.repeat(_branching/DP, DP)
        branch.append(weighted_data)
        #axes[2].bar(barindex + barwidth + dist, zimmstockmeyer_branch, barwidth,  edgecolor='white', color=colors[i])
        #axes[2].set_ylim(0, 2)

        index+=1

        barindex +=1

    violin_parts = axes[0].violinplot(branch,showmeans=True, showextrema=True)#
    violin_parts_2 = axes[1].violinplot(zimmstockmeyer_branching,showmeans=True, showextrema=True)#


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
    if 'cmeans' in violin_parts:
        violin_parts['cmeans'].set_color('black')
    # Customize colors

    for i, pc in enumerate(violin_parts_2['bodies']):
        pc.set_facecolor("grey")
        pc.set_edgecolor("grey")
        pc.set_alpha(0.7)


    # Customize range lines (extrema) if they exist
    if 'cmins' in violin_parts_2:
        violin_parts_2['cmins'].set_color('black')
    if 'cmaxes' in violin_parts:
        violin_parts_2['cmaxes'].set_color('black')
    if 'cbars' in violin_parts:
        violin_parts_2['cbars'].set_color('black')
    if 'cmeans' in violin_parts:
        violin_parts_2['cmeans'].set_color('black')


    ####### FOR THE BRANCHING FACTOR
    means_bf = [segment[0, 1] for segment in violin_parts['cmeans'].get_segments()]

    data_bf = {
        'Group': ['Grasses'] * 4 + ['Hardwoods'] * 7 + ['Softwoods'] * 5,
        'BranchingParameter': means_bf[:16]
    }

    df_bf = pd.DataFrame(data_bf)
# Perform ANOVA
    anova_result = stats.f_oneway(
        df_bf[df_bf['Group'] == 'Grasses']['BranchingParameter'],
        df_bf[df_bf['Group'] == 'Hardwoods']['BranchingParameter'],
        df_bf[df_bf['Group'] == 'Softwoods']['BranchingParameter']
    ) 
    print(f"F-statistic: {anova_result.statistic}, p-value: {anova_result.pvalue}")
    # Perform Tukey's HSD test
    tukey_result = pairwise_tukeyhsd(df_bf['BranchingParameter'], df_bf['Group'], alpha=0.05)
    print("Test results for branching factor")
    print(tukey_result) 

    ####### FOR THE G3
    means_g3 = [segment[0, 1] for segment in violin_parts_2['cmeans'].get_segments()]
    data_g3 = {
        'Group': ['Grasses'] * 4 + ['Hardwoods'] * 7 + ['Softwoods'] * 5,
        'BranchingParameter': means_g3[:16]
    }

    df_g3 = pd.DataFrame(data_g3)
# Perform ANOVA
    anova_result = stats.f_oneway(
        df_g3[df_g3['Group'] == 'Grasses']['BranchingParameter'],
        df_g3[df_g3['Group'] == 'Hardwoods']['BranchingParameter'],
        df_g3[df_g3['Group'] == 'Softwoods']['BranchingParameter']
    ) 
    print(f"F-statistic: {anova_result.statistic}, p-value: {anova_result.pvalue}")
    # Perform Tukey's HSD test
    tukey_result = pairwise_tukeyhsd(df_g3['BranchingParameter'], df_g3['Group'], alpha=0.05)
    print("Test results for g3")
    print(tukey_result) 
    

    # Define the legend labels and colors
    legend_labels = ["Experimental data", "C: All parameters free, limited S/G"]
    legend_colors = ["k", color3]

    # Create custom legend handles using Patch for stacked bars
    legend_handles = [patches.Patch(color=color, label=label) for color, label in zip(legend_colors, legend_labels)]

    #bond_distribution_MONADDSG, bond_distribution_ALLFREE, bond_distribution_ALLFREELIMSG,
    # Create a legend for the main plot
    #legend = fig.legend(handles=legend_handles, labels=legend_labels, loc='lower center', bbox_to_anchor=(0.6, -0.2), fontsize=fontsize, ncol=2)
    # Adjust the layout to prevent overlapping titles    
    axes[1].tick_params(axis='x', labelsize=fontsize)  # Change the fontsize of x-axis ticks
    axes[0].tick_params(axis='y', labelsize=fontsize)  # Change the fontsize of y-axis ticks    axes[0].tick_params(axis='x', labelsize=fontsize)  # Change the fontsize of x-axis ticks
    axes[1].tick_params(axis='y', labelsize=fontsize)  # Change the fontsize of y-axis ticks
        

    plt.tight_layout()
    savestring = "compare_branching.png"
    plt.savefig(savestring, bbox_inches='tight')
    return()


def calculate_alpha(pdi):
    pd_inv = 1/pdi
    if pd_inv < 0.63:
        alpha = (-0.58*pd_inv) + 0.36
    else:
        alpha = 0.05
    return(alpha)

def plot_PDI_branching(library_list, biomasses, labels):
    fontsize = 22
    fig, axes = plt.subplots(2,figsize=(16,7), sharex=True)
     # Set the axis labels and title
    axes[1].set_xlabel('Biomass', fontsize=fontsize)
    axes[1].set_ylabel('Branching factor $\\alpha$', fontsize=fontsize)
    axes[0].set_ylabel('PDI', fontsize=fontsize)

    # Bar position and shape properties 
    dist = 0.0
    subtr = -0.05
    separation = 0.2
    width = 0.65  # Reduce the width for the target distribution
    barwidth = 0.22
    cmap = plt.get_cmap('gnuplot2')
    colors = [cmap(i) for i in np.linspace(0, 1, 12)]
    color1 = colors[4]
    color2 = colors[6]
    color3 = colors[8]
    barindex=0
    index = 0
    labs =[]

    color4 = colors[2]
    colorlist= [color1, color2, color3]
    # Set positions of each biomass


    plotlablist = ["A", "B", "C"]

    #new_euclidian_distances =  np.zeros(len(new_libraries))
    #old_euclidian_distances =  np.zeros(len(old_libraries))

    N = 500 

    transparency = 1
    for i in range(len(library_list[0])):
        if 'MN' in biomasses[index] and 'MW' in biomasses[index]:
            labs.append(labels[index])
            print(i)
            print("All params free")
            mn_ALLFREE = calculate_number_avg(library_list[0][index][1][:N])
            mw_ALLFREE = calculate_weight_avg(library_list[0][index][1][:N])
            pd_ALLFREE = mw_ALLFREE/mn_ALLFREE
            alpha_ALLFREE =calculate_alpha(pd_ALLFREE)
            #pd_ALLFREE = np.average(DPs)#, weights = DPs)
            print("All params free limsg")

            mn_ALLFREELIMSG = calculate_number_avg( library_list[1][index][1][:N])
            mw_ALLFREELIMSG = calculate_weight_avg(library_list[1][index][1][:N])
            pd_ALLFREELIMSG = mw_ALLFREELIMSG/mn_ALLFREELIMSG
            alpha_ALLFREELIMSG = calculate_alpha(pd_ALLFREELIMSG)

            print("MonaddSG")
            mn_MONADDSG = calculate_number_avg(library_list[2][index][1][:N])
            mw_MONADDSG = calculate_weight_avg(library_list[2][index][1][:N])
            pd_MONADDSG = mw_MONADDSG/mn_MONADDSG
            alpha_MONADDSG = calculate_alpha(pd_MONADDSG)


            alpha_list = [alpha_MONADDSG, alpha_ALLFREE, alpha_ALLFREELIMSG]
            pd_list = [pd_MONADDSG, pd_ALLFREE, pd_ALLFREELIMSG]
            for sim in range(len(library_list)):

                axes[1].bar(barindex-1+ separation*sim + barwidth + dist, alpha_list[sim], barwidth,  edgecolor='white', color=colorlist[sim])

                #axes[2].text(barindex-1+ separation*sim + barwidth + dist, -0.25, plotlablist[sim], fontsize=12,horizontalalignment = 'center')
                axes[0].bar(barindex-1+ separation*sim + barwidth + dist, pd_list[sim], barwidth,  edgecolor='white', color=colorlist[sim])

            #if 'MN' in biomasses[index] and 'MW' in biomasses[index]:
                #target_DP = (biomasses[index]["MN"]*(biomasses[index]["Monolignols"]["S"]/100)/S_monolignol_weight) + (biomasses[index]["MN"]*(biomasses[index]["Monolignols"]["G"]/100)/G_monolignol_weight)
                target = biomasses[index]["MW"]/biomasses[index]["MN"]
                axes[0].bar(barindex-1 , target, barwidth - subtr,  edgecolor='white', color="k")

                target = calculate_alpha(target)
            #if 'MN' in biomasses[index] and 'MW' in biomasses[index]:

                if 'alpha' in biomasses[index]:
                    #target_DP = (biomasses[index]["MN"]*(biomasses[index]["Monolignols"]["S"]/100)/S_monolignol_weight) + (biomasses[index]["MN"]*(biomasses[index]["Monolignols"]["G"]/100)/G_monolignol_weight)
                    target = biomasses[index]["alpha"]#["MN"]/biomasses[index]["MW"]
                if target <0.05:
                    target = 0.01
                axes[1].bar(barindex -1, target, barwidth - subtr,  edgecolor='white', color="k")
            barindex +=1



        index +=1

    # Define the legend labels and colors
    legend_labels = ["Experimental data","A: Monomer addition rate and S/G ratio free", "B: All parameters free", "C: All parameters free, limited S/G"]
    legend_colors = ["k", color1, color2, color3]

    x_positions = np.arange(len(labs))
    axes[1].set_xticks(x_positions-1 + 0.25)
    axes[1].set_xticklabels(labs, fontsize=fontsize)

    # Create custom legend handles using Patch for stacked bars
    legend_handles = [patches.Patch(color=color, label=label) for color, label in zip(legend_colors, legend_labels)]



    #bond_distribution_MONADDSG, bond_distribution_ALLFREE, bond_distribution_ALLFREELIMSG,
    # Create a legend for the main plot
    legend = fig.legend(handles=legend_handles, labels=legend_labels, loc='lower center', bbox_to_anchor=(0.5, -0.2), fontsize=fontsize, ncol=2)
    # Adjust the layout to prevent overlapping titles    
    axes[0].tick_params(axis='x', labelsize=fontsize)  # Change the fontsize of x-axis ticks
    axes[0].tick_params(axis='y', labelsize=fontsize)  # Change the fontsize of y-axis ticks
    axes[1].tick_params(axis='y', labelsize=fontsize)  # Change the fontsize of y-axis ticks    axes[0].tick_params(axis='x', labelsize=fontsize)  # Change the fontsize of x-axis ticks
        

    plt.tight_layout()
    savestring = "compare_PDI_branching.png"
    plt.savefig(savestring, bbox_inches='tight')

def plot_parameter_correlation(library_list, fitting_parameters_list):
    """
    Plots the Pearson correlation matrix between input fitting parameters and output bond distribution data.

    Parameters:
    -----------
    library_list : list
        A nested list containing bond distribution data. Each entry in the list represents a different set of 
        libraries, with each library containing a list of bonds and their distribution.
    
    fitting_parameters_list : list
        A list containing fitting parameters. Each entry in the list is a sub-list representing the fitting parameters 
        corresponding to the libraries in `library_list`.

    Returns:
    --------
    None
        The function generates and saves a heatmap figure representing the Pearson correlation matrix. The heatmap 
        visually displays the correlation between input fitting parameters and output bond distributions, with color 
        indicating the strength and direction of the correlations.
    """
    fig, ax = plt.subplots(figsize=(12,10))
    fontsize = 22
    
    inputs = np.array(fitting_parameters_list[0])  #
    inputs = np.delete(inputs, [1, 2], axis=1)#np.vstack((fitting_parameters_list[0],fitting_parameters_list[1],fitting_parameters_list[2]))#

    outputs =[]
    for j in range(1):
        for i in range(len(library_list[0])):
            bond_distribution = calculate_avg_bond_distribution(library_list[j][i][1])
            bonds = bond_distribution.values()
            #print(bonds)
            outputs.append(list(bonds)[:])

    outputs = np.array(outputs)
    print(np.shape(outputs))
    print(np.shape(inputs))
    import pandas as pd
    # Combine the arrays into a single DataFrame
    data = np.hstack((inputs, outputs))
    columns = ["Mon add\nrate", "S/G ratio", "G-G", "S-S", "S-G", "$\\beta-O-4$", "$\\beta-\\beta$", "$\\beta-5$", "$5-5$", "$\\alpha-O-4$", "$4-O-5$", "$\\beta-1$"]
    #columns = ['A', 'B', 'C', 'W', 'X', 'Y', 'Z', 'A', 'B', 'C', 'W', 'X', 'Y', 'Z']
    df = pd.DataFrame(data, columns=columns)

    # Compute the correlation matrix
    correlation_matrix = df.corr(method='pearson')
        # Set the axis labels and title
    ax.set_xlabel('Input parameters     Output parameters', fontsize=fontsize)
    ax.set_ylabel('Input parameters     Output parameters', fontsize=fontsize)
    ax.tick_params(axis='x', labelsize=fontsize)  # Change the fontsize of x-axis ticks
    ax.tick_params(axis='y', labelsize=fontsize)  # Change the fontsize of y-axis ticks
    # Display the correlation matrix
    print(correlation_matrix)
    #correlation = stats.pearsonr(inputs,outputs)
    ticks = [-1, -0.5, 0, 0.5, 1]
    
    sns.heatmap(correlation_matrix, cmap='coolwarm_r', cbar_kws={'ticks': ticks, "label": "Pearson correlation"}, vmin=-1, vmax=1)
         # Set the axis labels and title    
    
    cbar = ax.collections[0].colorbar
    ax.collections[0].colorbar.set_ticklabels([t for t in ticks], fontsize=fontsize)
    cbar.set_label('Pearson correlation', fontsize=fontsize)

    ax.axhline(5, color="k", linestyle="--", linewidth=2)
    ax.axvline(5, color="k", linestyle="--", linewidth=2)

    savestring = "compare_corr.png"
    # Adjust the layout to prevent overlapping titles    
    #plt.ylim(len(columns),0)
# Reverse the y-axis
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(savestring, bbox_inches='tight')


# Function to get the sorting key
def get_sort_key(dict_name):
    return biomasses[dict_name]["short label"]
biomasses = load_biomass_data("biomasses.json")
sorted_dict_names = sorted(biomasses.keys(), key=get_sort_key)

location_string = "example_data/set_C/"

target_data = []
monaddsglibraries = []
allparamslibs = []
allparamslimitedsglibs =[]

fitting_params_list = []
fitting_data = []
for biomass_label in list(sorted_dict_names):
    library_import_string = location_string + biomass_label + "/BEST_FIT/best_Run/Output/library.json"
    if os.path.exists(library_import_string):
        allparamslimitedsglibs.append(import_library(library_import_string, "ligninkmc_"))
        fitting_import_string = location_string + biomass_label + "/best_kin_specs.txt"
        fitting_data.append(np.loadtxt(fitting_import_string))
        target_data.append(biomasses[biomass_label])
    else: print(library_import_string)

fitting_params_list.append(fitting_data)
fitting_data = []

location_string = "example_data/set_B/"
for biomass_label in list(sorted_dict_names):
    library_import_string = location_string + biomass_label + "/BEST_FIT/best_Run/Output/library.json"
    if os.path.exists(library_import_string):
        allparamslibs.append(import_library(library_import_string, "ligninkmc_"))
        fitting_import_string = location_string + biomass_label + "/best_kin_specs.txt"
        fitting_data.append(np.loadtxt(fitting_import_string))

    else: print(library_import_string)

fitting_params_list.append(fitting_data)
fitting_data = []

location_string = "example_data/set_A/"
for biomass_label in list(sorted_dict_names):
    library_import_string = location_string + biomass_label + "/BEST_FIT/best_Run/Output/library.json"
    if os.path.exists(library_import_string):
        monaddsglibraries.append(import_library(library_import_string, "ligninkmc_"))
        fitting_import_string = location_string + biomass_label + "/best_kin_specs.txt"
        fitting_data.append(np.loadtxt(fitting_import_string))

    else: print(library_import_string)
fitting_params_list.append(fitting_data)

labels = [dict["short label"] for dict in target_data]

lib = [allparamslibs, allparamslimitedsglibs, monaddsglibraries]
plot_PDI_branching(lib, target_data, labels)
#plot_compare_SG_ratios_barchart(lib, target_data, labels)
#plot_g3_branchingfactor_subplots(lib, target_data, labels)
#plot_branching_subplots(lib, target_data, labels)
#plot_polydispersityindex_comparison(lib, target_data, labels)
#plot_compare_SGratio_freeenergies__heatmap(fitting_params_list, target_data, labels)
#plot_compare_bond_distributions(lib, target_data, labels)
#plot_parameter_correlation(lib, fitting_params_list)
#plot_compare_monomeradditionrates(lib, target_data, fitting_params_list, labels)
#plot_compare_SGratio_freeenergies(fitting_params_list, target_data, labels)
#plot_compare_SGratio_freeenergies_diff(fitting_params_list, target_data, labels)
#plot_mass_comparison(lib, target_data, labels)
#plot_energybarrier_change_barcharts(fitting_params_list, labels)
## plot_energybarrier_change(fitting_params_list, labels)
#plot_compare_euclidian_distances(lib, target_data, labels)
##plot_compare_SG_ratios(allparamslibs, allparamslimitedsglibs, monaddsglibraries, target_data, labels)
