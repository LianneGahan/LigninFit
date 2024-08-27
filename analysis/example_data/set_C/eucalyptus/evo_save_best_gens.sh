#!/bin/bash

#$1 : number of generations


rm -r BEST_FIT 2>/dev/null

mkdir BEST_FIT

touch BEST_FIT/var_vs_gen.txt



i=1
while [[ $i -le $1 ]]
do

	echo "$i 	 "$(sort -k 2,2 family_1/vars_"$i".txt | head -n 1)"" >> BEST_FIT/var_vs_gen.txt
	((i = i + 1))

done


touch BEST_FIT/lowest_vars_gen.txt

sort -k 3,3 BEST_FIT/var_vs_gen.txt | head -1 >> BEST_FIT/lowest_vars_gen.txt
echo "BEST FIT"

echo "$(sort -k 3,3 BEST_FIT/var_vs_gen.txt | head -1)"

echo "Copying best Gen"


while IFS=" " read gen_no run_no gen_var
do
	cp -r family_1/Generation_"$gen_no"/best_Run/ BEST_FIT
    
done < "BEST_FIT/lowest_vars_gen.txt"
