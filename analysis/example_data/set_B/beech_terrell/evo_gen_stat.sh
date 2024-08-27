#!/bin/bash
#$1: input parameter of the number of generations in the fitting algorithm
#$2: input parameter 

path="family_1/vars_*.txt"
num_files=$(ls -1 $path 2>/dev/null | wc -l)

#echo echo -n "" > statoos.txt

if [[ $num_files -eq $1 ]]
then
	echo "START_GEN			0" >> generations.log
else
	to_del=`expr $num_files + 1`
	echo "LAST_BROKEN_GEN			$to_del" >> generations.log
	echo "removing interrupted Generation number: $to_del"
	rm -r family_1/Generation_$to_del
fi
