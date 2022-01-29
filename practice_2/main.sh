#!/bin/bash

############################################################
############################################################
# Main program                                             #
############################################################
############################################################

PATH_DATA="/home/data/practica2_genomica/"

for file in "${PATH_DATA}*.fastq.gz"
do
    echo "Running MetaPhlAn on ${file}"

    filename=$(basename ${file})
    BOWTIE_RESULT="bowtie2/${filename}.bowtie2.bz2"
    PROFILE_RESULT="profiles/${filename}_profiled.tsv"

    # Launch Metaphlan
    metaphlan ${file} --input_type fastq --bowtie2out ${BOWTIE_RESULT} -o ${PROFILE_RESULT}
done