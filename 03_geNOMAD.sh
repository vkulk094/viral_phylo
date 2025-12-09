#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=0                   # give the job all node memory
#SBATCH --time=11:59:00
#SBATCH --ntasks-per-node=48
#SBATCH --job-name=geNOMAD
#SBATCH --output=geNOMAD_%A_%a.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email_here@domain.com

## OUTPUT ------------------------------------------------------------------
## Results written under: $SCRATCH/genomad_results/<biome>_<region>/<sample>_genomad
## -------------------------------------------------------------------------
# Strict bash: fail on errors, undefined vars, and pipeline failures
set -euo pipefail
IFS=$'\n\t'

## USAGE -------------------------------------------------------------------
## sbatch --array=0-180%181 geNOMAD.sh
## Requires: contig_list.txt with full paths to final.contigs.fa
## -------------------------------------------------------------------------

# Parameterize list file (override with: sbatch -e LISTFILE=other.txt ...)
LISTFILE="${LISTFILE:-contig_list.txt}"
[[ -f "$LISTFILE" ]] || { echo "ERROR: List file '$LISTFILE' not found" >&2; exit 1; }

# Fetch the FASTA path for this array task
FASTA_IN=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$LISTFILE")
[[ -s "$FASTA_IN" ]] || { echo "ERROR: Missing or empty FASTA at index $SLURM_ARRAY_TASK_ID" >&2; exit 1; }

# Derive sample, region, biome from path hierarchy
sample=$(basename "$(dirname "$FASTA_IN")")
region=$(basename "$(dirname "$(dirname "$FASTA_IN")")")
biome=$(basename "$(dirname "$(dirname "$(dirname "$FASTA_IN")")")")

# User-defined variables
GENOMAD_DB="/home/vaibhavk/scratch/genomad_db"
OUTDIR="$SCRATCH/genomad_results/${biome}_${region}/${sample}_genomad"
######################################

# Load modules dependencies
module load StdEnv/2023 cudacore/.12.6.2 aragorn/1.2.41 mmseqs2/17-b804f python/3.10

# Activate your virtual environment in $SLURM_TMPDIR
. ~/geNOMAD_env/bin/activate

# Create output directory and log start
mkdir -p "$OUTDIR"
echo "[$(date)] Processing $FASTA_IN → $OUTDIR"

# Run geNOMAD end-to-end pipeline with cleanup
genomad end-to-end \
        --cleanup \
        "$FASTA_IN" \
        "$OUTDIR" \
        "$GENOMAD_DB"
    
