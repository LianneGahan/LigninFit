 

## Instructions for how to use the paramater optimization algorithm


### Setup directory and data to be fitted

- Create a directory in which you want to run the simulations (this will be referred to as topdir)

- copy the content of this directory into topdir

- If you have changed the simulation code, copy the new code into topdir/latest/, and run compile.sh via the terminal command

./compile.sh

Check that the "Params" directory, and the "Output" directory exist within the directory "latest"

- create an empty directory called "family_1" in topdir.

- Open the file "best_kin_specs.txt" and enter those kinetic parameters from which the optimization algorithm should start. This file will be updated continuously during the algorithm

- IMPORTANT: Within topdir, there is a file "keywords.txt" (if not, you need to create it). From this file the algorithm takes
the names of the experimental data sets, and therefore it needs to contain AT LEAST ONE keyword. This is now automated and the code writes into the "keywords.txt" file, reading from the experimental data contained in "latest/Output/expe_data/expe_saccharification_keyword_glc.txt" or "latest/Output/expe_data/expe_saccharification_keyword_xyl.txt". If you want to put in the keywords manually you need to comment in evo_all_in_one.sh: line 16: ./evo_keywords.sh.


- Open the file "kin_params_to_randomize.txt". Here, you can choose, which parameters should be fitted by setting the respective value to 0 (fixed at provided value) or 1 (fitted, provided value taken as starting point). 

- Open the files "max_kin_vals.txt" and "max_grid_search.txt". Here, you can specify the maximum value for each parameter.

- Open the files "min_kin_vals.txt" and "min_init_vals.txt". Here, you can specify the minimum value for each parameter.


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

- Now, to start the fitting algorithm simply execute the shell script "evo_wrapper.sh" via `./evo_wrapper.sh`.

- It also keeps track of the completed generations and any interruptions to the fitting algorithm due to server restart, broken pipe (if not run in background) etc in the file 'generations.log'.

- To resume from an interruption, simply run `./evo_wrapper.sh` and it will pick up again from the interrupted generation.

- While the algorithm runs, the files "best_candidate_specs.txt" are updated if the difference between simulated and experimental data is reduced compared to the previous set of parameters in the file. Furthermore, the current lowest error is contained in the file "smallest_var.txt"

- The best fitted set of parameters are saved in the directory BEST_FIT/best_run/Params. The resulting lignin library generated under these conditions is saved in "BEST_FIT/best_Run/Output/library.json".

- If you want to visually compare the simulated and experimental data, run the code in a separate folder and use the parameters found in in the directory BEST_FIT/best_run/Params (see Code README at the top of this repository for more detail).

- The 'goodness' of the fit is determined by taking the Euclidian distance value comparing the model output average bond distribution and experimental data average bond distribution. 
