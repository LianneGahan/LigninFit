#!/bin/bash

#$1 : number of families
#$2 : number of generations
#$3 : number of folders in each generation
#$4 : DELTA
#$5 : Max number of cores to be used
#$6 : Resume from generation
#$7 : Number of Samples to use in testing
#(echo "99999999">smallest_var.txt)



(echo "0">no_reduction_count.txt)
#(echo "Generation 1")
(python3 reset_gradient.py)



(for i in $(seq $6 $2)
do

    (echo "Generation $i")
	./evo_create_folders.sh $1 ${i} $3
	PROCESSID=$!
	wait $PROCESSID
	(python3 rando_specs.py $i $3 $4 $1)
	PROCESSID=$!
	wait $PROCESSID
    (echo "Running")
	./evo_run.sh $1 $i $3 $5 $7
	PROCESSID=$!
	wait $PROCESSID
    (echo "Averaging")
	./evo_average_fit_and_print_variance.sh $1 $i $3 $5
	PROCESSID=$!
	wait $PROCESSID
	rm vars.txt
	./evo_pick_min.sh $1 $i $3 
	PROCESSID=$!
	wait $PROCESSID
    python3 find_min_vars.py $1 $i
	PROCESSID=$!
	wait $PROCESSID
	cp vars.txt family_$1/vars_$i.txt
	PROCESSID=$!
	wait $PROCESSID
	./evo_clear_gen.sh $1 $i $(<best_candidate.txt)
	PROCESSID=$!
	wait $PROCESSID
	echo "LAST_COMPLETED_GEN		$i" >> generations.log

done
)
echo "ALL_GOOD_NOW			0" >> generations.log

echo "Choosing best Gen"

./evo_save_best_gens.sh $2
