# Constants declaration
SRC="./src"
SCRIPT_BOWTIE_ALIGN="${SRC}/script_align_bowtie.sh"

# Load scripts
source $SCRIPT_BOWTIE_ALIGN

# Constants declaration
REFERENCE_GENOME="s_cerevisiae.fasta"
NUMBER_MISSMATCHES=1

# Launch script allowing 1 missmatch
main_bowtie_align $REFERENCE_GENOME $NUMBER_MISSMATCHES