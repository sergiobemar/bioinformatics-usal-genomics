#!/bin/bash

############################################################
############################################################
# Main program                                             #
############################################################
############################################################

# Constants declaration
SRC="./src"
SCRIPT_COLLECT_DATA="${SRC}/script_collect_data.sh"
SCRIPT_BOWTIE_ALIGN="${SRC}/script_align.sh"
SCRIPT_UTILS="${SRC}/utils.sh"
PATH_DATA="./data"
PATH_RESULTS="./results"
PATH_BOWTIE_INDEX="${PATH_DATA}/bowtie_index"

# Arguments for alignment
REFERENCE_GENOME="${PATH_DATA}/s_cerevisiae.fasta"

# Remove folders
rm -rf $PATH_BOWTIE_INDEX

# Create folders
mkdir -p  $PATH_DATA
mkdir -p $PATH_RESULTS

# Load scripts
source $SCRIPT_UTILS
#source $SCRIPT_COLLECT_DATA
source $SCRIPT_BOWTIE_ALIGN

# Create index from the reference gen
# Directory for index
mkdir $PATH_BOWTIE_INDEX

# Create index files variable and launch the index generation
BOWTIE_INDEX_FILES="${PATH_BOWTIE_INDEX}/index_`basename ${REFERENCE_GENOME%.*}`"

function_bowtie_create_index $REFERENCE_GENOME $BOWTIE_INDEX_FILES

# Launch script allowing 1 missmatch
NUMBER_MISMATCHES=1
echo "$(function_get_now) launching command: main_bowtie_align ${BOWTIE_INDEX_FILES} ${NUMBER_MISMATCHES} ${PATH_DATA} ${PATH_RESULTS}"
main_bowtie_align $BOWTIE_INDEX_FILES $NUMBER_MISMATCHES $PATH_DATA $PATH_RESULTS

# Launch script allowing 3 missmatch
NUMBER_MISMATCHES=3
echo "$(function_get_now) launching command: main_bowtie_align ${BOWTIE_INDEX_FILES} ${NUMBER_MISMATCHES} ${PATH_DATA} ${PATH_RESULTS}"
main_bowtie_align $BOWTIE_INDEX_FILES $NUMBER_MISMATCHES $PATH_DATA $PATH_RESULTS