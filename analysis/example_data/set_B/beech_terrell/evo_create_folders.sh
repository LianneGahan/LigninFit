#!/bin/bash

#$1 : family
#$2 : Number of this generation
#$3 : number of folders in each generation
#for i in $(seq 1 $2)
#do
mkdir family_$1/Generation_$2
for j in $(seq 1 $3)
do
	mkdir family_$1/Generation_$2/Run_${j}
	cp -a latest/Code family_$1/Generation_$2/Run_${j}/
	cp -a latest/Output/ family_$1/Generation_$2/Run_${j}/
	cp -a latest/Params/ family_$1/Generation_$2/Run_${j}/
done
#done

