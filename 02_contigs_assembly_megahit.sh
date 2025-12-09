#!/bin/bash

# Path to the file containing SRA IDs to be processed
SRA_IDS_FILE="sra_ids.txt"
INITIAL_DIR=$(pwd)
CURRENT_DIR_NAME=$(basename "$PWD")

# Check if the SRA IDs file exists
if [[ ! -f $SRA_IDS_FILE ]]; then
    echo "Error: File $SRA_IDS_FILE not found!"
    exit 1
fi

# Function to create a full job script for an SRA ID (single-node full job)
create_full_job_script() {
    local SRA_ID="$1"
    local JOB_SCRIPT="job_${SRA_ID}.sh"
    cat <<EOT > $JOB_SCRIPT
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=0
#SBATCH --time=11:59:59
#SBATCH --job-name=${CURRENT_DIR_NAME}_${SRA_ID}
#SBATCH --output=${CURRENT_DIR_NAME}_${SRA_ID}_%j.log

# Load necessary modules
module load StdEnv/2020
module load gcc/9.3.0
module load sra-toolkit/3.0.0
module load trimmomatic/0.39
module load megahit/1.2.9

# Set SRA_ID
SRA_ID="${SRA_ID}"

# Capture the initial working directory
INITIAL_DIR=\$(pwd)
OUTPUT_DIR="\${INITIAL_DIR}/${SRA_ID}_output"
TRIMMOMATIC_ADAPTER="/home/vaibhavk/adapters/TruSeq3-PE.fa"
LOCKFILE="\${INITIAL_DIR}/${SRA_ID}.lock"

# Skip if lock file exists
if [ -e "\$LOCKFILE" ]; then
    echo "Skipping \$SRA_ID as it is already being processed."
    exit 0
fi

# Create lock file
touch "\$LOCKFILE"

# Create necessary directories
mkdir -p \$OUTPUT_DIR

# Step 1: Download and convert SRA to FASTQ
cd \$OUTPUT_DIR
echo "Downloading and converting SRA to FASTQ for \$SRA_ID..."
prefetch --max-size 100GB \$SRA_ID
if [ \$? -ne 0 ]; then
    echo "Error: prefetch failed for \$SRA_ID"
    rm -f "\$LOCKFILE"
    exit 1
fi
fasterq-dump --split-files \$SRA_ID -O .
if [ \$? -ne 0 ]; then
    echo "Error: fasterq-dump failed for \$SRA_ID"
    rm -f "\$LOCKFILE"
    exit 1
fi

# Step 2: Trimming with Trimmomatic
echo "Trimming FASTQ files for \$SRA_ID..."
java -jar \$EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -basein \${SRA_ID}_1.fastq -baseout \${SRA_ID}.fastq \
    ILLUMINACLIP:\$TRIMMOMATIC_ADAPTER:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36
if [ \$? -ne 0 ]; then
    echo "Error: Trimmomatic failed for \$SRA_ID"
    rm -f "\$LOCKFILE"
    exit 1
fi

# Step 3: Assemble with MEGAHIT using 48 threads (single node)
echo "Assembling trimmed reads with MEGAHIT for \$SRA_ID..."
if [ -d "megahit_\${SRA_ID}" ]; then
    echo "Resuming MEGAHIT assembly for \$SRA_ID..."
    megahit --continue -t 48 -o megahit_\${SRA_ID}
else
    megahit -1 \${SRA_ID}_1P.fastq -2 \${SRA_ID}_2P.fastq -t 48 -o megahit_\${SRA_ID}
fi
if [ \$? -ne 0 ]; then
    echo "Error: MEGAHIT failed for \$SRA_ID"
    rm -f "\$LOCKFILE"
    exit 1
fi

# Step 4: Cleanup intermediate files after success
echo "Cleaning up intermediate files for \$SRA_ID..."
rm -f \${SRA_ID}_1.fastq \${SRA_ID}_2.fastq \${SRA_ID}.sra

echo "Pipeline completed successfully for \$SRA_ID."
rm -f "\$LOCKFILE"
cd \$INITIAL_DIR
EOT
}

# Function to create a MEGAHIT resume job script for an SRA ID (two-node resume job)
create_resume_job_script() {
    local SRA_ID="$1"
    local JOB_SCRIPT="job_${SRA_ID}.sh"
    cat <<EOT > $JOB_SCRIPT
#!/bin/bash
#SBATCH --job-name=${CURRENT_DIR_NAME}_MEG_${SRA_ID}
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=48
#SBATCH --mem=0
#SBATCH --time=23:59:59
#SBATCH --output=MEG_${CURRENT_DIR_NAME}_${SRA_ID}_%j.log

# Load necessary modules
module load StdEnv/2020
module load gcc/9.3.0
module load megahit/1.2.9

# Set SRA_ID
SRA_ID="${SRA_ID}"

# Capture the initial working directory
INITIAL_DIR=\$(pwd)
OUTPUT_DIR="\${INITIAL_DIR}/${SRA_ID}_output"
LOCKFILE="\${INITIAL_DIR}/${SRA_ID}.lock"

# Path to the checkpoint directory
CHECKPOINT_DIR="\${OUTPUT_DIR}/megahit_${SRA_ID}/"

cd \${OUTPUT_DIR}
pwd

# Resume the MEGAHIT job using 96 threads (2 nodes × 48 cores)
megahit --continue -t 96 -o \${CHECKPOINT_DIR}
MEGAHIT_EXIT_CODE=\$?

if [ \$MEGAHIT_EXIT_CODE -ne 0 ]; then
    echo "Error: MEGAHIT resume failed for \$SRA_ID"
    rm "\${LOCKFILE}"
    exit 1
fi

echo "MEGAHIT completed successfully for \$SRA_ID."
rm "\${LOCKFILE}"
cd \$INITIAL_DIR
EOT
}

# Loop through each SRA_ID from the file (processing in reverse order)
tac "$SRA_IDS_FILE" | while IFS= read -r SRA_ID; do
    # Define the lock file and checkpoint directory paths
    LOCKFILE="${INITIAL_DIR}/${SRA_ID}.lock"
    CHECKPOINT_DIR="${INITIAL_DIR}/${SRA_ID}_output/megahit_${SRA_ID}/"
    
    # Check if final contigs or a done file exist in the checkpoint directory
    if [[ -e "${CHECKPOINT_DIR}/final.contigs.fa" || -e "${CHECKPOINT_DIR}/done" ]]; then
        echo "Skipping $SRA_ID as it is already completed."
        continue
    fi

    # Skip if lock file exists
    if [[ -e "$LOCKFILE" ]]; then
        echo "Skipping $SRA_ID as it is already being processed."
        continue
    fi

    # Decide which job script to create: resume if checkpoints exist, else full job
    if [[ -e "${CHECKPOINT_DIR}/checkpoints.txt" ]]; then
        echo "Checkpoint found for $SRA_ID. Creating resume job script."
        create_resume_job_script "$SRA_ID"
    else
        echo "No checkpoint for $SRA_ID. Creating full job script."
        create_full_job_script "$SRA_ID"
    fi

    # Submit the generated SLURM script
    sbatch job_${SRA_ID}.sh
    # Remove the temporary SLURM script after submission
    rm job_${SRA_ID}.sh
done


