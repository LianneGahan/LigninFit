#!/bin/bash

# $1 = family
# $2 = Generation number
# $3 = Number to save

echo "Cleaning house";
cd family_$1
mv Generation_$2/Run_$3 Run_$3
rm -r Generation_$2/*
mv Run_$3 Generation_$2/best_Run
touch Generation_$2/best_was_run_number_$3.txt
cd ../

echo "Done cleaning house";


