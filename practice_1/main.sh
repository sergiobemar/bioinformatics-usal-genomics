# Constants declaration
SRC="./src"
SCRIPT_COLLECT_DATA="${SRC}/script_collect_data.sh"
SCRIPT_BOWTIE_ALIGN="${SRC}/script_align_bowtie.sh"
PATH_DATA="./data"
PATH_RESULTS="./results"
PATH_BOWTIE_INDEX="${PATH_DATA}/bowtie_index"

# Arguments for alignment
REFERENCE_GENOME="${PATH_DATA}/s_cerevisiae.fasta"

# Load scripts
#source $SCRIPT_COLLECT_DATA
source $SCRIPT_BOWTIE_ALIGN

# Create index from the reference gen
main_bowtie_create_index $REFERENCE_GENOME $PATH_BOWTIE_INDEX

# Launch script allowing 1 missmatch
NUMBER_MISSMATCHES=1
main_bowtie_align $PATH_BOWTIE_INDEX $NUMBER_MISSMATCHES $PATH_DATA $PATH_RESULTS
