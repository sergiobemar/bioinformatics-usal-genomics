#!/bin/bash

##############################################################################
# 
# Main program which taking fastq.gz files as input, launch metaphlan tool
# in order to get the composition of microbial communities from the different.
#
# In addition, BowTie files are saved too so that it will be faster to relaunch
# if it was needed.
#
# Arguments:
#   - Directory where the fastq.gz files are stored in
#
# More info in the following links:
#	- GitHub repository MetaPhlAn 3.0: https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-3.0
#	- MetaPhlAn website: https://huttenhower.sph.harvard.edu/metaphlan
#
##############################################################################

PATH_DATA=$1
PATH_FASTQ_FILES="${PATH_DATA}/*.fastq.gz"

PATH_BOWTIE_FOLDER="bowtie2"
PATH_PROFILES_FOLDER="profiles"

# Create folders if not exists
if [[ ! -e $PATH_BOWTIE_FOLDER ]]; then
	echo "$(date +'%Y-%m-%d %H:%M:%S'): Creating folder ${PATH_BOWTIE_FOLDER}"
	mkdir -p $PATH_BOWTIE_FOLDER
fi


if [[ ! -e $PATH_PROFILES_FOLDER ]]; then
	echo "$(date +'%Y-%m-%d %H:%M:%S'): Creating folder ${PATH_PROFILES_FOLDER}"
	mkdir -p $PATH_PROFILES_FOLDER
fi

for file in $PATH_FASTQ_FILES
do

    echo ""
    echo "$(date +'%Y-%m-%d %H:%M:%S'): Running MetaPhlAn on ${file}"

    # Get files' names
    filename=$(basename ${file})
    BOWTIE_RESULT="${PATH_BOWTIE_FOLDER}/${filename}.bowtie2.bz2"
    PROFILE_RESULT="${PATH_PROFILES_FOLDER}/${filename}_profiled.tsv"

    # Launch Metaphlan command
    echo "$(date +'%Y-%m-%d %H:%M:%S'): Launching command: metaphlan ${file} --input_type fastq --bowtie2out ${BOWTIE_RESULT} -o ${PROFILE_RESULT}"
    metaphlan ${file} --input_type fastq --bowtie2out ${BOWTIE_RESULT} -o ${PROFILE_RESULT}
done

echo "$(date +'%Y-%m-%d %H:%M:%S'): Execution finished"
