#!/bin/bash

# Remove contents of keywords.txt file
echo -n "" > keywords.txt

#directory="latest/Output/expe_data"

#if ls "$directory/expe_saccharification_"*_glc.txt >/dev/null 2>&1; then
    # Loop through matching file names for glc files
#    for glc_file in latest/Output/expe_data/expe_saccharification_*_glc.txt; do
#        # Extract wild card entry from glc file name
#        entry=$(echo $glc_file | awk -F'_' '{print $4}')
#        # Write entry to keywords.txt file if it doesn't already exist
#        if ! grep -Fxq "$entry" keywords.txt; then
#            echo $entry >> keywords.txt
#        fi
#    done

#fi

#if ls "$directory/expe_saccharification_"*_xyl.txt >/dev/null 2>&1; then
#    # Loop through matching file names for xyl files
#    for xyl_file in latest/Output/expe_data/expe_saccharification_*_xyl.txt; do
#        # Extract wild card entry from xyl file name
#        entry=$(echo $xyl_file | awk -F'_' '{print $4}')
#        # Write entry to keywords.txt file if it doesn't already exist
#        if ! grep -Fxq "$entry" keywords.txt; then
#            echo $entry >> keywords.txt
#        fi
#    done
#fi
