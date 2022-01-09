# Constants declaration
SRC="./src"
SCRIPT_COLLECT_DATA="${SRC}/script_collect_data.sh"
SCRIPT_BOWTIE_ALIGN="${SRC}/script_align_bowtie.sh"
PATH_DATA="./data"
PATH_RESULTS="./results"
PATH_BOWTIE_INDEX="${PATH_DATA}/bowtie_index"

# Arguments for alignment
REFERENCE_GENOME="${PATH_DATA}/s_cerevisiae.fasta"

# Remove folders
rm -rf $PATH_BOWTIE_INDEX

# Load scripts
#source $SCRIPT_COLLECT_DATA
source $SCRIPT_BOWTIE_ALIGN

# Create index from the reference gen
# Directory for index
mkdir $PATH_BOWTIE_INDEX

# Create index files variable and launch the index generation
INDEX_FILES="${PATH_BOWTIE_INDEX}/index_`basename ${REFERENCE_GENOME%.*}`"

function_bowtie_create_index $REFERENCE_GENOME $INDEX_FILES

# Launch script allowing 1 mismatch
NUMBER_MISMATCHES=1
main_bowtie_align $PATH_BOWTIE_INDEX $NUMBER_MISMATCHES $PATH_DATA $PATH_RESULTS
