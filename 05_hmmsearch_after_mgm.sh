#!/bin/bash

# Directories and tools
BASE_DIR="/home/vaibhavk/scratch/final_contigs/DNA_Folder"  #Change to desired DNA folder
MGM_DIR="/home/vaibhavk/MetaGeneMark_linux_64/mgm"
VIRAL_PROT="/home/vaibhavk/scratch/viral_prot"
OUTPUT_DIR="$BASE_DIR/MGM_Output"   # Final output directory inside DNA_Folder
INTRM_DIR="$BASE_DIR/MGM_INTRM"     # Intermediate files directory inside DNA_Folder
SLURM_LOG_DIR="$INTRM_DIR/slurm_logs"
TIME="2:55:00"          # Wall time (2 hours and 55 minutes)

# Ensure directories exist
mkdir -p "$OUTPUT_DIR" "$INTRM_DIR" "$SLURM_LOG_DIR"

# Submit combined MetaGeneMark and HMMER jobs
echo "[Optimized] Submitting combined MGM and HMMER jobs..."
find "$BASE_DIR" -type d -name 'megahit_SR*' | while read -r folder; do
    ORIGINAL_NAME=$(basename "$folder")
    SRA_ID=${ORIGINAL_NAME#megahit_}
    CONTIGS_FILE="$folder/final.contigs.fa"
    MGM_OUTPUT="$INTRM_DIR/DNA_contigs_mgm_${SRA_ID}.faa"
    MGM_LIST="$INTRM_DIR/sequence_${SRA_ID}.lst"
    HMM_OUTPUT="$OUTPUT_DIR/${SRA_ID}_hmm.out"

    # Skip if contigs file does not exist
    if [ ! -f "$CONTIGS_FILE" ]; then
        echo "Skipping $SRA_ID: no final.contigs.fa found."
        continue
    fi

    # Write job script for combined MetaGeneMark and HMMER
    JOB_FILE="$INTRM_DIR/combined_${SRA_ID}.sbatch"
    cat <<- EOF > "$JOB_FILE"
#!/bin/bash
#SBATCH --job-name=mgm_hmmer_$SRA_ID
#SBATCH --output=$SLURM_LOG_DIR/combined_${SRA_ID}.out
#SBATCH --error=$SLURM_LOG_DIR/combined_${SRA_ID}.err
#SBATCH --time=$TIME
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=0

# Load required modules
module load hmmer/3.4

# Step 1: Run MetaGeneMark
echo "[SLURM] Running MetaGeneMark for $SRA_ID..."
"$MGM_DIR/gmhmmp" -m "$MGM_DIR/MetaGeneMark_v1.mod" -A "$MGM_OUTPUT" -o "$MGM_LIST" "$CONTIGS_FILE"

# Check if MetaGeneMark succeeded
if [ ! -s "$MGM_OUTPUT" ]; then
    echo "Error: MetaGeneMark output is empty for $SRA_ID. Exiting job."
    exit 1
fi

# Step 2: Run HMMER
echo "[SLURM] Running HMMER for $SRA_ID..."
hmmsearch -E 1e-10 --cpu=$NTASKS "$VIRAL_PROT/clusters/vir_models.hmm" "$MGM_OUTPUT" > "$HMM_OUTPUT"

# Check if HMMER succeeded
if [ ! -s "$HMM_OUTPUT" ]; then
    echo "Error: HMMER output is empty for $SRA_ID. Exiting job."
    exit 1
fi

echo "[SLURM] Finished processing for $SRA_ID."
EOF

    # Submit the job to SLURM
    sbatch "$JOB_FILE"
done

echo "[Checkpoint] All jobs submitted to SLURM."
