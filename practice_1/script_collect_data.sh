#!/bin/bash

##############################################################################
# Return the current date tiem in string format without spaces, for namefiles
##############################################################################
function_get_now_logs() {
    date +"%Y%m%d_%H%M%S"
}

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

copy_sample_files() {


    # Copy original gen samples (fastq.gz) to my folder 
    for f in $(find /home/practicasGenomica/ -name "*.fastq.gz")
    do
        # Get all names, removing extension and convert to lower
        file=$(basename $f)
        file=$(echo "${file%.fastq.gz}")
        folder=$(echo $file |  awk '{print tolower($0)}')

        # Create folder for fastaq files
        echo ""
        echo "--------------------------------------------------------------------------------------------------------------------------------------"
        echo "$(function_get_now) copy_sample_files(): Creating directory: ${folder}"
        echo "--------------------------------------------------------------------------------------------------------------------------------------"
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

index_ref_gen() {

    # Copy the reference genome to our path
    cp /home/practicasGenomica/S_cerevisiae.fasta ./

    # Index the reference genome using bowtie-build
    echo ""
    echo "--------------------------------------------------------------------------------------------------------------------------------------"
    echo "$(function_get_now) index_ref_gen(): Start indexing to ${1}"
    echo "--------------------------------------------------------------------------------------------------------------------------------------"
    bowtie-build S_cerevisiae.fasta $1
}

main() {
    # Copy the sample files
    copy_sample_files

    # Create the index
    index_ref_gen "index_ref"
}

main