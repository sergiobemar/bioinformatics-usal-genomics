#!/bin/bash


##############################################################################
# Return the current date tiem in string format without spaces, for namefiles
##############################################################################
function_get_now() {
    date +"%Y-%m-%d %H:%M:%S"
}

##############################################################################
# Rename all the files of the input directory using lowercase
#
# Arguments:
#   - Input directory
##############################################################################
rename_to_lower_files() {

    echo ""
    echo "--------------------------------------------------------------------------------------------------------------------------------------"
    echo "$(function_get_now) rename_to_lower_files(): Renaming to lower the files in directory ${1}"
    echo "--------------------------------------------------------------------------------------------------------------------------------------"
    for i in $( ls ${1} | grep [A-Z] )
    do
        mv -i $i `echo $i | tr 'A-Z' 'a-z'`
    done
}

##############################################################################
# Copy all sample of genomes to the currect directory.
# All the files are in the directory /home/practicasGenomica.
#
# In addition to copy the files, this functions extracts all the copied fastq.gz
##############################################################################
copy_sample_files() {


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
        folder=$(echo $file |  awk '{print tolower($0)}')

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
    gunzip -r ./
}

##############################################################################
# Copy reference genome to the current directory
#
# Arguments:
#   - Path of the reference genome
##############################################################################
copy_reference_genome() {
    
    # Input parameters to more readable variables
    input_reference_genome=$1

    # Create the output index genome reference
    output_file=$(echo "./`basename $input_reference_genome`" | tr 'A-Z' 'a-z')

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
#
# Arguments:
#   - Path of the reference genome
##############################################################################
main() {

    # Input parameters to more readable variables
    input_reference_genome=$1

    # Copy the sample files
    copy_sample_files

    # Copy the reference_genome
    copy_reference_genome $input_reference_genome
}

# Constants declaration
REFERENCE_GENOME_PATH="/home/practicasGenomica/S_cerevisiae.fasta"

# Call main method
main $REFERENCE_GENOME_PATH
