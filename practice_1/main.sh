# Constants declaration
SRC="./src"
SCRIPT_COLLECT_DATA="${SRC}/script_collect_data.sh"
SCRIPT_BOWTIE_ALIGN="${SRC}/script_align_bowtie.sh"
DATA="./data"

# Arguments for alignment
REFERENCE_GENOME="s_cerevisiae.fasta"

# Load scripts
#source $SCRIPT_COLLECT_DATA
source $SCRIPT_BOWTIE_ALIGN

# Launch script allowing 1 missmatch
NUMBER_MISSMATCHES=1
main_bowtie_align $REFERENCE_GENOME $NUMBER_MISSMATCHES $DATA