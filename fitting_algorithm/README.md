## Instructions for how to use the paramater optimization algorithm
The figure shows the scheme of the computational framework LigninFit. In the leftmost box, the "Initial conditions identification" finds the starting parameters for the "Fitting procedure" by comparing the target experimental bond distribution to an existing simulated database of bond distributions, following a grid search. In a loop, the "Fitting procedure" optimises the input parameters (S/G ratio, monomer addition rate, and the free energy barriers of S-S, G-G, and S-G bond formation) such that the simulated bond distribution matches the target experimental one.
Following the identification of the best fit, a large library of molecules is generated for the corresponding set of parameter values, leading to the "Compilation of molecule information" and the identification of the "Optimal growth conditions" (rightmost box).
![image info](./overview.png)

The sections below describe the setup and use of LigninFit. The main code can be found under "fitting_algrithm", and the companion analysis codes to produce the figures from the publication (in preparation) can be found in the "analysis" directory. 

**Please Note:** The fitting algorithm requires Lignin-KMC to function correctly, so please ensure you download a copy of the codes from the relevant GitHub directory and include it as specified below.

### Setup directory and data to be fitted

- Create a directory in which you want to run the simulations (this will be referred to as topdir)

- Copy the content of this directory into topdir

- Download Lignin-KMC and copy "lignin-kmc" repository into topdir/latest/Code/


Check that the "Params" directory, and the "Output" directory exist within the directory "latest"

- create an empty directory called "family_1" in topdir.

- Open the file "best_kin_specs.txt" and enter those kinetic parameters from which the optimization algorithm should start. This file will be updated continuously during the algorithm

- Open the file "kin_params_to_randomize.txt". Here, you can choose, which parameters should be fitted by setting the respective value to 0 (fixed at provided value) or 1 (fitted, provided value taken as starting point). 

- Open the files "max_kin_vals.txt" and "max_grid_search.txt". Here, you can specify the maximum value for each parameter.

- Open the files "min_kin_vals.txt". Here, you can specify the minimum value for each parameter.


### Running the algorithm

- Once the setup is complete, you may run the algorithm. For this, you need to specify the following five parameters:

    - The number of families **N_family** (currently to be kept at a value of 1)
    - The number of generations per family **N_gens** (choose an integer greater than or equal to 1)
    - The number of subsets per generation **N_subsets** (choose an integer greater than or equal to 1)
    - The percentage around which to vary the parameters **delta** (choose a number between 0 and 1)
    - The number of CPU cores you have available **N_cores** (choosing a number higher than the actual number of cores available will reduce the efficiency of the algorithm considerably)
    - The Number of repeats in each subset **N_repeats**
    - The number of samples per degree of freedom for the initial gridsearch **N_samples** (Keep in mind that this will scale M**N_samples, where M is the number of degrees of freedom in the gridsearch) 
    
- The above seven values or fit settings are taken from the file 'fit_settings.txt'. The values are tab separated entries in the file. A default set of values are provided. But can be changed depending on resources available and the 'closeness' of the initial starting paramter set.

- You need to also set the 3 parameters in "topdir/Params/simulation_parameters.txt"

    - The maximum simulation time (seconds) **max_time** (set to a very high value e.g. 100000000 to ensure polymerisation is more complete)
    - A bool to determine additional levels of filewriting **filewriting**
    - A numeric gridsearch option **gridsearch_option**
        - **0** Ranges are provided for the **N** parameters that are chosen to be varied, as well as the number of different values (noted **N_samples**) they each can take during this step of Initial conditions identification. Then, Lignin-KMC simulations are performed for these $N^M$ sets of parameter values, and the closest matching bond distribution (with the lowest Euclidean distance) that respects the range of values provided for each parameter is identified. Ranges are provided for the parameters that are chosen to be varied. These can be the S/G ratio (i.e. the ratio of S to G monolignols that are added to the simulation system), the monomer addition rate (i.e. the rate at which monomers are added to the simulation system), the S-S, G-G, and S-G  scaling parameters (i.e. each parameter scales specific free energy barriers for bond formation, which are converted into interaction rates within Lignin-KMC), and the minimum and maximum number of monomers (i.e. the number of monolignols in the system when the simulation starts and the maximum it can contain, respectively) A new gridsearch is performed using the limits "min_kin_vals.txt" and "gridsearch_max.txt". 
        - **1** to use pre-existing gridsearch data "topdir/gridsearch/gridsearch.json". The bond distribution of the target experimental biomass is compared to pre-existing simulation results stored in a database. Within this database, the closest matching bond distribution (with the lowest Euclidean distance) that respects the range of values provided for each parameter ("min_kin_vals.txt" and "max_kin_vals.txt") is identified.
        - **2** No search is performed, and the initial conditions for the Fitting procedure (in the middle box) are manually entered. The user chooses the starting point established in "topdir/Params/kinetic_params.txt"

- Now, to start the fitting algorithm simply execute the shell script "evo_wrapper.sh" via `./evo_wrapper.sh`.

- It also keeps track of the completed generations and any interruptions to the fitting algorithm due to server restart, broken pipe (if not run in background) etc in the file 'generations.log'.

- To resume from an interruption, simply run `./evo_wrapper.sh` and it will pick up again from the interrupted generation.

- While the algorithm runs, the files "best_candidate_specs.txt" are updated if the difference between simulated and experimental data is reduced compared to the previous set of parameters in the file. Furthermore, the current lowest error is contained in the file "smallest_var.txt"

- The best fitted set of parameters are saved in the directory BEST_FIT/best_run/Params. The resulting lignin library generated under these conditions is saved in "BEST_FIT/best_Run/Output/library.json".

- If you want to visually compare the simulated and experimental data, run the code in a separate folder and use the parameters found in in the directory BEST_FIT/best_run/Params (see Code README at the top of this repository for more detail).

- The 'goodness' of the fit is determined by taking the Euclidian distance value comparing the model output average bond distribution and experimental data average bond distribution. 
