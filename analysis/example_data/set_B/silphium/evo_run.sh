#!/bin/bash

# $1 = family
# $2 = Generation number
# $3 = Number of runs in generation
# $4 = Maximum number of cores available
# $5 = Number of simulations in Sample 


COUNT=1
PIDCOUNT=1
echo "Running code";
for i in $(seq 1 $3)
do
	cd family_$1/Generation_$2/Run_${COUNT}
	python3 Code/simulation_main.py $5 8 #>/dev/null &
	sleep 3 
	pids[${PIDCOUNT}]=$!	
	PIDCOUNT=$((${PIDCOUNT} + 1));
		
	joblist=($(jobs -p))

	while [ ${#joblist[*]} -ge $4 ]
	do
		sleep 0.1
		joblist=($(jobs -p))
	done

	cd ../../../
	let COUNT++
done

for pid in ${pids[*]}; do
    wait $pid
done
echo "Done running code";
wait $PROCESSID


