#!/bin/bash

# Load libraries
source ./utils.sh

##############################################################################
# Alignes a genome with its referency gen using bowtie 
# with the number of missmatches setting by an argument.
#
# It creates a .sam file with the alignement which will be compressed to a .bam file
#
# Arguments:
#   - Reference gen indexed previously
#   - Number of missmatches
#   - Path of directory where the data is stored
##############################################################################
function_bowtie_alignement_genome() {

    # Use a more readable variables for input parameters
    index_ref=$1
    number_missmatch=$2
    path_data=$3

    # Create directory for results
    results_directory="${path_data}/bowtie_results_missmatch_missmatches_${number_missmatch}"
    mkdir $results_directory

    for f in $(find $path_data -name "*.fastq")
    do
       
        # Extract file (with extension) and filename (without extension)
        file=$(basename $f)
        filename_result="$results_directory/${file%.*}_${2}"
        dir=$(dirname $f)

        # Create the variables with the names for files SAM, BAM and BAI
        sam_result="${filename_result}.sam"
        bam_result="${filename_result}.bam"
        bai_result="${filename_result}.bai"

        echo ""
        echo "--------------------------------------------------------------------------------------------------------------------------------------"
        echo "$(function_get_now) function_alignement_genome(): Align file ${f} to ${sam_result}"
        echo "--------------------------------------------------------------------------------------------------------------------------------------"

        # Alignement using bowtie with 1 missmatch
        echo ""
        echo "$(function_get_now) function_alignement_genome(): Aligning genome using bowtie"
        #bowtie -t -p 16 -v $number_missmatch -S $index_ref $f $sam_result

        # Compress .sam to .bam file
        #function_compress_sam $sam_result $bam_result

        # Index .bam result to baiecho ""
        echo "$(function_get_now) function_alignement_genome(): Indexing BAM file, the result will be saved in ${bai_result}"
        #samtools index $bam_result $bai_result
    done
}

##############################################################################
# Compress a .sam file input to a .bam file. Then, it remove the .sam file
#
# Both arguments must have their extension
#
# Arguments:
#   - Input .sam file
#   - Name of the output .bam file
##############################################################################
function_compress_sam() {
    echo ""
    echo "-------------------------------------------------------------------"
    echo "$(function_get_now) function_compress_sam(): Compressing SAM file ${1} to BAM file ${2}"
    echo "-------------------------------------------------------------------"
    samtools sort -@ 16 -o $2 $1
    rm $1
}

##############################################################################
# Index a genome in fasta format using bowtie library
#
# Arguments:
#   - Name of the output index gen
#   - Input genome where is the fasta file
##############################################################################
function_bowtie_create_index() {

    # Use a more readable variables for input parameters
    output_index_name=$1
    input_genome_path=$2

    # Index the reference genome using bowtie-build
    echo ""
    echo "--------------------------------------------------------------------------------------------------------------------------------------"
    echo "$(function_get_now) function_bowtie_create_index(): Start indexing ${input_genome_path} to ${output_index_name}"
    echo "--------------------------------------------------------------------------------------------------------------------------------------"
    bowtie-build $input_genome_path $output_index_name
}

##############################################################################
# Main function to lauch all the script to align the sample genes with
# the index allowing n missmatches passed as argument
#
# Arguments:
#   - Path of the index gen
#   - Number of missmatches
##############################################################################
main_bowtie_align() {

    # Use a more readable variables for input parameters
    input_reference_genome=$1
    number_missmatch=$2

    # Create directory for index
    path_index="bowtie_index"
    mkdir $path_index

    # Index reference genome
    path_data="data"
    index_files="${path_data}/${path_index}/index_`basename ${input_reference_genome%.*}`"
    function_bowtie_create_index $index_files $input_reference_genome

    # Align gen samples with reference genome
    function_bowtie_alignement_genome $index_files $number_missmatch $path_data
}
