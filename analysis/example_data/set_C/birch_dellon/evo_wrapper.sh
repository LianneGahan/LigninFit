#!/bin/bash

#param1 : number of families (currently fixed to 1)
#param2 : number of Gens in family
#param3 : number of subfolders in a Gen
#param4 : DELTA
#param5 : number of cores to be used
#param6 : Resume interrupted fitting from this generation number
#param7 : Number of runs to use in a sample
#param8 : CHOICE: 0 to carry out a new gridsearch, 1 to use pre-existing gridsearch data, 
#			2 for the user to choose the starting point
#param9 : Number of samples per degree of freedom in the gridsearch** ignored if param8 == 0

# read 5 tab-separated input parameters from fit_params.txt
read -r param1 param2 param3 param4 param5 param7 param9 < fit_settings.txt
read -r paramA paramB paramC < latest/Params/simulation_parameters.txt
# read a interruption status from generations.log
#read -r param6 < statoos.txt
last_gen=$(tail -1 "generations.log" | awk '{print $2}')

param6=`expr $last_gen + 1`
#param6=$(cut -f2 generations.log)

# print the values of the input parameters
echo "--------------------------------------------"
echo "number of families: $param1"
echo "Generations per family: $param2"
echo "Sub-folders per generation: $param3"
echo "DELTA: $param4"
echo "cores to be used: $param5"
echo "Number of runs to use in a sample: $param7"
echo "Grid search option: $paramC"
echo "Grid search degrees of freedom: $param9 (ignored if Grid search option is not 0)"
echo "--------------------------------------------"


if [[ $last_gen -gt 0 ]]
then
	echo "Resume fitting from generation: $param6"
	./evo_gen_stat.sh $param2
	PROCESSID=$!
	wait $PROCESSID
	./evo_all_in_one.sh $param1 $param2 $param3 $param4 $param5 $param6 $param7
else
	./evo_janitor.sh
	PROCESSID=$!
	wait $PROCESSID
	./evo_find_starting_conditions.sh $param5 $param7 $paramC $param9
	PROCESSID=$!
	wait $PROCESSID
	./evo_all_in_one.sh $param1 $param2 $param3 $param4 $param5 1 $param7
fi

./evo_janitor.sh


# Do the final run with the found kinetic parameters
cd BEST_FIT/best_Run/
echo "Running final big simulation!"
python3 Code/simulation_main.py 200 $param5 >/dev/null &
# Calculate final costs
cd ../../
#python3 evaluate_best_fit.py

#echo -e '\033[1m\033[39mALL DONE!! \033[0m'
echo -e 'Find results in the folder >>\033[1m\033[35m BEST_FIT \033[0m'
