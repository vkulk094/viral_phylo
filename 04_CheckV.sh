#!/bin/bash

#SBATCH --ntasks=12
#SBATCH --mem-per-cpu=2G
#SBATCH --time=01:15:00       
#SBATCH --job-name=CheckV-terres
#SBATCH --output=logs/CheckV_%A_%a.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email_here@vaibhavk.com

## USAGE -------------------------------------------------------------------
## sbatch --array=1-158 CheckV.sh
## Requires: genomad_virus_list.txt with full paths to <sample>_genomad/final.contigs_summary/final.contigs_virus.fna
## -------------------------------------------------------------------------
## OUTPUT ------------------------------------------------------------------
## Results written under: $SCRATCH/checkv_results/<biome>_<region>/<sample>_checkv
## -------------------------------------------------------------------------

# Strict bash: fail on errors, undefined vars, and pipeline failures
set -euo pipefail
IFS=$'\n\t'

##########  USER VARIABLES  ###############################################
# Parameterize list file (override with: sbatch -e LISTFILE=other.txt ...)
LISTFILE="${LISTFILE:-genomad_virus_list.txt}"
[[ -f "$LISTFILE" ]] || { echo "ERROR: List file '$LISTFILE' not found" >&2; exit 1; }

mkdir -p logs

# Determine viral FASTA either from array list or positional argument
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    VIRAL_FASTA=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LISTFILE")
elif [[ -n "${1:-}" ]]; then
    VIRAL_FASTA="$1"
else
    echo "ERROR: No viral FASTA specified."
    echo "  Standalone usage: sbatch CheckV.sh /path/to/<sample>_virus.fna"
    exit 1
fi
[[ -s "$VIRAL_FASTA" ]] || { echo "ERROR: FASTA '$VIRAL_FASTA' missing or empty" >&2; exit 1; }

# Derive sample, region, biome from path hierarchy
sample=$(basename "$(dirname "$(dirname "$VIRAL_FASTA")")" _genomad)
base=$(basename "$(dirname "$(dirname "$(dirname "$VIRAL_FASTA")")")")  # e.g. marine_eq
biome=${base%%_*}
region=${base##*_}

CHECKV_DB="/home/vaibhavk/scratch/checkv-db-v1.5"
OUTDIR="$SCRATCH/checkv_results/${biome}_${region}/${sample}_checkv"
THREADS=${SLURM_CPUS_ON_NODE:-1}
###########################################################################

# Load modules dependencies
module load StdEnv/2020 hmmer/3.3.2 prodigal-gv/2.6.3 diamond/2.0.4 python/3.10

# Generate your virtual environment in $SLURM_TMPDIR
. ~/CheckV_env/bin/activate

# Create output directory and log start
mkdir -p "$OUTDIR"
echo "[$(date)] Running CheckV on $VIRAL_FASTA → $OUTDIR"

# Run CheckV end-to-end pipeline
checkv end_to_end \
        "$VIRAL_FASTA" \
        "$OUTDIR" \
        -t "$THREADS" \
        -d "$CHECKV_DB"
    