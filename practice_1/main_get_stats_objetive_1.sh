#!/bin/bash

############################################################
############################################################
# Main program                                             #
############################################################
############################################################

# Constants declaration
SRC="./src"
SCRIPT_UTILS="${SRC}/utils.sh"

# Params
PATH_RESULTS=$1
RESULTS_FILE=$2

# Load scripts
source $SCRIPT_UTILS

# Create file for results
RESULTS_FILE="${PATH_RESULTS}/$RESULTS_FILE"

# Iterate over .bam files to launch the stats
for f in $(find $PATH_RESULTS -name "*.bam")
do
    echo ""
    echo "$(function_get_now): Stats for BAM file: ${f}"
    
    echo "" >> $RESULTS_FILE
    echo "--------------------------------------------------------------------------------------------------------------------------------------" >> $RESULTS_FILE
    echo "$(function_get_now): Stats for BAM file: ${f}" >> $RESULTS_FILE
    bamtools stats -in $f >> $RESULTS_FILE
    echo "--------------------------------------------------------------------------------------------------------------------------------------" >> $RESULTS_FILE
done