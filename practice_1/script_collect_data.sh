#!/bin/bash

rename_to_lower_files() {

    echo "rename_to_lower_files(): Renaming to lower the files in directory ${1}"
    for i in $( ls ${1} | grep [A-Z] )
    do
        mv -i $i `echo $i | tr 'A-Z' 'a-z'`
    done
}

copy_sample_files() {


    # Copy original gen samples (fastq.gz) to my folder 
    for f in $(find /home/practicasGenomica/ -name "*.fastq.gz")
    do
        file=$(basename $f)
        file=$(echo "${file%.fastq.gz}")
        folder=$(echo $file |  awk '{print tolower($0)}')

        # Create folder for fastaq files
        echo "copy_sample_files(): Creating directory: ${folder}"
        mkdir $folder

        # Copy files to directory
        echo "copy_sample_files(): Copy file ${f} to ${folder}"
        cp $f $folder

        rename_to_lower_files $(find $folder -name "*.fastq.gz")
    done

    # Uncompress all fastq.gz files
    echo "copy_sample_files(): Uncompressing files"
    gunzip -r ./
}

index_ref_gen() {

    cp /home/practicasGenomica/S_cerevisiae.fasta ./

    echo "index_ref_gen(): Start indexing to ${1}"
    bowtie-build S_cerevisiae.fasta $1
}

copy_sample_files

index_ref_gen "index_ref"
