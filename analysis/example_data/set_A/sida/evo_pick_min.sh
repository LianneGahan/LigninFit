#!/bin/bash

# $1 : family
# $2 : Generation
# $3 : Runs in generation
COUNT=$2
MIN=0
touch vars.txt
cd family_$1/Generation_${COUNT}/Run_1/Output/
MIN=$(<var.txt)||MIN=999999999999
echo 1 $MIN >> ../../../../vars.txt
cd ../../../../
for i in $(seq 2 $3)
do
	cd family_$1/Generation_${COUNT}/Run_${i}/Output/
	MIN=$(<var.txt) || MIN=999999999999
	#MIN=$(<(( $MIN <= $(<var.txt)) && echo "$a" || echo "$b"))
	echo $MIN
	echo $i		$MIN >> ../../../../vars.txt
	cd ../../../../
done


