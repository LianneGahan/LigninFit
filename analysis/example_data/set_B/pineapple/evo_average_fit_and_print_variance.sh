#!/bin/bash

# $1 : Family
# $2 : Generation
# $3 : Runs in generation
# $4 : Max number of CPUs
COUNT=$2
FAM=$1
GEN=$2
RUNS=$3
#RUNS=$((${RUNS}+1))
echo "Runs = ${RUNS}"
PIDCOUNT=1
#readarray -t names < keywords.txt
#for i in $(seq 1 $#)
#do
#	names[$i] = $1
#	echo "$1 and ${${names}[$i]}"
#	shift
#done

for i in $(seq 1 ${RUNS})
do
	echo "in evo_average_fit_and_print_variance: i = $i"
    python3 average_fit_and_print_variance.py family_${FAM}/Generation_${COUNT}/Run_${i}/ #>/dev/null &
	pids[${PIDCOUNT}]=$!
	PIDCOUNT=$((${PIDCOUNT} + 1));	
	joblist=($(jobs -p))
	while [ ${#joblist[*]} -ge $4  ]
	do
		sleep 0.1
		joblist=($(jobs -p))
	done
done

for pid in ${pids[*]}; do
    wait $pid
done
wait $PROCESSID


