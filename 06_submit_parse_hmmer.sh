#!/bin/bash
#SBATCH --job-name=hmmer_all
#SBATCH --output=slurm_logs/hmmer_%A_%a.out
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48    # Set to match the number of CPU cores per node
#SBATCH --mem=0                   # give the job all node memory
#SBATCH --tmp=120G
#SBATCH --array=1-121

set -euo pipefail
mkdir -p slurm_logs

########################
#  Modules and env
########################
module load StdEnv/2020 r/4.2.2 python/3.10.2 mafft/7.471 hmmer/3.3.2 trimal/1.4 fasttree/2.1.11

# Make threads visible to tools that look for CPUS_PER_TASK
export SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-$SLURM_NTASKS}
THREADS=${SLURM_CPUS_PER_TASK}
export OMP_NUM_THREADS=$THREADS
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

# Keep tool scratch on node-local
export TMPDIR="${SLURM_TMPDIR:-/tmp}"
export MAFFT_TMPDIR="$TMPDIR"

echo "Node: $(hostname)"
echo "NTASKS: $SLURM_NTASKS  THREADS: $THREADS"
echo "SLURM_TMPDIR: ${SLURM_TMPDIR:-<unset>}"
df -hT "${SLURM_TMPDIR:-/tmp}" || true

########################
#  Manifest line
########################
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /home/vaibhavk/scratch/all_samples.txt)
read -r region sra_id _ <<< "$line"
echo "Region: $region  SRA: $sra_id"

########################
#  Paths (persistent)
########################
BASE_IN="/home/vaibhavk/scratch/checkv_results/${region}"
HMM_OUT="${BASE_IN}/MGM_Output/${sra_id}_hmm.out"
MGM_FAA="${BASE_IN}/MGM_INTRM/${sra_id}.faa"
BIGFAA="/home/vaibhavk/scratch/IMG_VR/all_clusters.faa"   # or .faa.gz

for f in "$HMM_OUT" "$MGM_FAA" "$BIGFAA"; do
  [[ -s "$f" ]] || { echo "Missing input: $f" >&2; exit 2; }
done

########################
#  Stage-in to local scratch
########################
WORK="$SLURM_TMPDIR/hmmer_${sra_id}"
mkdir -p "$WORK"
ulimit -n 4096 || true

# Copy per-SRA inputs
cp -f "$HMM_OUT" "$WORK/"; sync
cp -f "$MGM_FAA" "$WORK/"; sync

# Cache BIGFAA once per node so multiple array tasks do not all copy it
# Detect the node-local mount from SLURM_TMPDIR (e.g., /localscratch on Cedar, /local on SHARCNET)
NODE_LOCAL_ROOT="$(dirname "$SLURM_TMPDIR")"
NODE_CACHE="${NODE_LOCAL_ROOT}/.imgvr_cache"

# If we cannot create a node-wide cache (no permission on mount root), use a per-job cache
if ! mkdir -p "$NODE_CACHE" 2>/dev/null; then
  echo "No write permission on $NODE_LOCAL_ROOT; using job-local cache."
  NODE_CACHE="$SLURM_TMPDIR/.imgvr_cache"
  mkdir -p "$NODE_CACHE"
fi

BIGFAA_LOCAL="$NODE_CACHE/$(basename "$BIGFAA")"

# Optional: show where we are caching
ls -ld "$NODE_LOCAL_ROOT" "$NODE_CACHE" || true

# Use flock if available; otherwise fall back to a simple copy
if command -v flock >/dev/null 2>&1; then
  LOCK="$NODE_CACHE/.bigfaa.lock"
  exec 9> "$LOCK" || true
  flock -n 9 || true
  if [[ ! -s "$BIGFAA_LOCAL" ]]; then
    echo "Caching $(basename "$BIGFAA") to $BIGFAA_LOCAL"
    cp -f "$BIGFAA" "$BIGFAA_LOCAL".partial && mv -f "$BIGFAA_LOCAL".partial "$BIGFAA_LOCAL"
  fi
  flock -u 9 || true
else
  [[ -s "$BIGFAA_LOCAL" ]] || cp -f "$BIGFAA" "$BIGFAA_LOCAL"
fi

# Symlink into this job's workspace
ln -sf "$BIGFAA_LOCAL" "$WORK/$(basename "$BIGFAA")"

########################
#  Stage-out helper + traps
########################
PERSIST_BASE="/home/vaibhavk/scratch/hmmer_align/${region}/${sra_id}"
stage_out() {
  echo "Staging out to $PERSIST_BASE"
  mkdir -p "$PERSIST_BASE"
  shopt -s nullglob
  mv "$WORK"/hits_*.tar.gz "$PERSIST_BASE"/ 2>/dev/null || true
  if [[ -d "$WORK/hits_trees" ]]; then
    tar -czf "$PERSIST_BASE/hits_trees.tar.gz" -C "$WORK" hits_trees || true
  fi
  [[ -f "$WORK/hmmer_processing.log" ]] && cp -f "$WORK/hmmer_processing.log" "$PERSIST_BASE"/ || true
  [[ -f "$WORK/${sra_id}_hmmer_out.RData" ]] && cp -f "$WORK/${sra_id}_hmmer_out.RData" "$PERSIST_BASE"/ || true
}
trap 'echo "Caught signal, staging out..."; stage_out; exit 143' TERM INT USR1

########################
#  Run (single step that uses all 8 CPUs)
########################
SCRIPT="/home/vaibhavk/scratch/06_parse_hmmer.R"

# IMPORTANT: -n 1 = one task; -c $SLURM_NTASKS = give that task all 8 CPUs
# Do NOT use --cpu-bind=cores here unless you also pass -c, or you will pin to 1 core.
srun -n 1 -c "$SLURM_NTASKS" Rscript "$SCRIPT" "$region" "$sra_id" "$WORK/$(basename "$BIGFAA")"

# Normal completion stage-out
stage_out
echo "Done."
