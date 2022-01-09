#!/bin/bash


# Load libraries
source ./utils.sh

##############################################################################
# Copy all sample of genomes to the currect directory.
# All the files are in the directory /home/practicasGenomica.
#
# In addition to copy the files, this functions extracts all the copied fastq.gz
#
# Arguments:
#   - Parent output directory
##############################################################################
copy_sample_files() {

    # Input parameters to more readable variables
    data_directory=$1

    # Copy original gen samples (fastq.gz) to my folder 
    for f in $(find /home/practicasGenomica/ -name "*.fastq.gz")
    do
        echo ""
        echo "--------------------------------------------------------------------------------------------------------------------------------------"
        echo "$(function_get_now) copy_sample_files(): File ${f}"
        echo "--------------------------------------------------------------------------------------------------------------------------------------"

        # Get all names, removing extension and convert to lower
        file=$(basename $f)
        file=$(echo "${file%.fastq.gz}")
        folder=$(echo "$data_directory/$file" |  awk '{print tolower($0)}')

        # Create folder for fastaq files
        echo "$(function_get_now) copy_sample_files(): Creating directory: ${folder}"
        mkdir $folder

        # Copy files to directory
        echo "copy_sample_files(): Copy file ${f} to ${folder}"
        cp $f $folder

        rename_to_lower_files $(find $folder -name "*.fastq.gz")
    done

    # Uncompress all fastq.gz files
    echo ""
    echo "--------------------------------------------------------------------------------------------------------------------------------------"
    echo "$(function_get_now) copy_sample_files(): Uncompressing files"
    echo "--------------------------------------------------------------------------------------------------------------------------------------"
    gunzip -r $data_directory
}

##############################################################################
# Copy reference genome to the current directory
#
# Arguments:
#   - Path of the reference genome
#   - Parent output directory
##############################################################################
copy_reference_genome() {
    
    # Input parameters to more readable variables
    input_reference_genome=$1
    data_directory=$2

    # Create the output index genome reference
    output_file=$(echo "./$data_directory/`basename $input_reference_genome`" | tr 'A-Z' 'a-z')

    # Copy the reference genome to our path
    echo ""
    echo "--------------------------------------------------------------------------------------------------------------------------------------"
    echo "$(function_get_now) copy_reference_genome(): Copying reference genome from ${input_reference_genome} to ${output_file}"
    echo "--------------------------------------------------------------------------------------------------------------------------------------"

    cp $input_reference_genome $output_file
}

##############################################################################
# Main function which make the following steps:
#   1. Copy all samples of the genes
#   2. Index the reference genome
##############################################################################
main_collect_data() {

    # Input parameters to more readable variables
    input_reference_genome="/home/practicasGenomica/S_cerevisiae.fasta"

    # Create directory for data
    data_directory="data"
    mkdir $data_directory

    # Copy the sample files
    copy_sample_files $data_directory

    # Copy the reference_genome
    copy_reference_genome $input_reference_genome $data_directory
}
