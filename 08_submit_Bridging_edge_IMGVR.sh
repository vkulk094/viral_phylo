#!/usr/bin/env bash
#SBATCH --job-name=br_patristic_csv
#SBATCH --array=0-5                 # 6 tasks, one per region
#SBATCH --nodes=1
#SBATCH --ntasks=48                 # use cluster default behavior
#SBATCH --mem=0
#SBATCH --time=36:00:00
#SBATCH --output=%x_%A_%a.out

set -euo pipefail
module load StdEnv/2020 r/4.2.2 || true

# --- Config ---
BASE="/home/vaibhavk/scratch/hmmer_align"
REGIONS=(marine_south marine_north marine_eq terrestrial_south terrestrial_north terrestrial_eq)
CLADE_RSCRIPT="/home/vaibhavk/scratch/scripts/08_Bridging_edge_IMGVR.R"

# --- Knobs ---
UNTAR_CHUNK="${UNTAR_CHUNK:-2000}"   # members per chunk
COPY_LIMIT="${COPY_LIMIT:-8}"        # concurrent gz copies per node
MIN_LINES="${MIN_LINES:-2}"

# --- Env hygiene ---
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export SHELL=/bin/bash
export TMPDIR="${SLURM_TMPDIR:-/tmp}"
mkdir -p "$TMPDIR"

# --- Binaries ---
RSCRIPT_BIN="$(command -v Rscript || true)"
if [[ -z "$RSCRIPT_BIN" ]]; then
  echo "FATAL: Rscript not found."
  exit 127
fi

# --- Region setup ---
REGION="${REGIONS[$SLURM_ARRAY_TASK_ID]}"
ROOT="$BASE/$REGION"

echo "Region: $REGION"
echo "Node:   $(hostname)"
echo "TMPDIR: ${TMPDIR}"

# Fan-out uses SLURM_NTASKS (cluster default param you’re using)
J="${SLURM_NTASKS:-}"
if [[ -z "$J" || "$J" -lt 1 ]]; then J=$(nproc || echo 1); fi

FAIL_LOG="$HOME/clade_failures_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt"
: > "$FAIL_LOG"

# ---------- helpers ----------
_throttle() {
  local kind="$1" limit="$2"
  local lockdir="${TMPDIR}/.clade_tokens_${kind}"
  mkdir -p "$lockdir"
  while true; do
    local n
    n=$(ls -1 "$lockdir" 2>/dev/null | wc -l | tr -d ' ')
    if (( n < limit )); then
      local tok="$lockdir/$$.$RANDOM.token"
      : > "$tok"
      echo "$tok"
      return 0
    fi
    sleep 0.2
  done
}
_release_token() {
  local token="$1"
  [[ -n "${token:-}" && -f "$token" ]] && rm -f "$token" || true
}

append_csv_no_dup_header() {
  # $1 src, $2 dest
  if [[ ! -s "$2" ]]; then
    cat "$1" >> "$2"
  else
    awk 'NR>1{print}' "$1" >> "$2"
  fi
}

# ---------- per-SRA worker ----------
process_one() {
  local tar_gz="$1"
  local parent sra final_out dir_local local_gz local_tar manifest chunk_prefix tmp_out

  parent="$(dirname "$tar_gz")"
  sra="$(basename "$parent")"
  final_out="${parent}/hits_trees_clade_distance_patristic.csv"

  # Resume
  local lines=0
  if [[ -f "$final_out" ]]; then
    lines=$(wc -l < "$final_out" 2>/dev/null || echo 0)
  fi
  if [[ -z "${FORCE:-}" && -s "$final_out" && "$lines" -ge "$MIN_LINES" ]]; then
    echo "[SKIP ] CLADE  $final_out (lines=$lines)"
    return 0
  fi

  # Local staging
  dir_local="${TMPDIR}/clade_${sra}_$$"
  mkdir -p "$dir_local"
  local_gz="${dir_local}/hits_trees.tar.gz"
  local_tar="${dir_local}/hits_trees.tar"
  manifest="${dir_local}/members.exact"
  chunk_prefix="${dir_local}/members.part."
  tmp_out="${dir_local}/hits_trees_clade_distance_patristic.csv"
  : > "$tmp_out"

  echo "[STAGE] Copy -> local (gz)  $tar_gz -> $local_gz"
  local tok_copy
  tok_copy=$(_throttle copy "$COPY_LIMIT")
  if ! cp -f "$tar_gz" "$local_gz"; then
    echo "[FAIL copy] $tar_gz" | tee -a "$FAIL_LOG"
    _release_token "$tok_copy"
    rm -rf "$dir_local" || true
    return 0
  fi
  _release_token "$tok_copy"

  echo "[STAGE] Gunzip to .tar      $local_gz -> $local_tar"
  if ! gzip -dc "$local_gz" > "$local_tar"; then
    echo "[FAIL gunzip] $tar_gz" | tee -a "$FAIL_LOG"
    rm -rf "$dir_local" || true
    return 0
  fi

  echo "[INFO ] Listing members (exact) -> $manifest"
  # Keep tar names verbatim from tar -tf. Only strip CR. Filter by extensions case-insensitively.
  local tmp_list="${dir_local}/members.raw"
  if ! tar -tf "$local_tar" > "$tmp_list"; then
    echo "[FAIL list] $tar_gz" | tee -a "$FAIL_LOG"
    rm -rf "$dir_local" || true
    return 0
  fi

  awk 'BEGIN{IGNORECASE=1}
       {gsub(/\r$/,"")}
       tolower($0) ~ /\.(tre|tree|trees|treefile|nwk|newick)$/ { print $0 }' \
       "$tmp_list" > "$manifest"

  # Split exact-name list into chunks
  if [[ "${UNTAR_CHUNK}" -gt 0 ]]; then
    split -l "${UNTAR_CHUNK}" -d -a 4 "$manifest" "$chunk_prefix"
  else
    cp -f "$manifest" "${chunk_prefix}0000"
  fi

  echo "[RUN  ] CLADE  $tar_gz -> stream chunks from $local_tar (chunk_size=${UNTAR_CHUNK})"

  # Iterate chunks and call R once per chunk, streaming members directly
  local part part_idx=0
  shopt -s nullglob
  for part in ${chunk_prefix}*; do
    local chunk_csv="${dir_local}/chunk_${part_idx}.csv"
    MEMBERS_FILE="$part" "$RSCRIPT_BIN" "$CLADE_RSCRIPT" "$local_tar" "$chunk_csv" || {
      echo "[FAIL clade] $tar_gz (chunk ${part_idx})" | tee -a "$FAIL_LOG"
      rm -f "$chunk_csv" || true
      part_idx=$((part_idx+1))
      continue
    }

    if [[ -s "$chunk_csv" ]]; then
      append_csv_no_dup_header "$chunk_csv" "$tmp_out"
    fi
    rm -f "$chunk_csv" || true
    part_idx=$((part_idx+1))
  done
  shopt -u nullglob

  # Validate and promote final CSV
  local new_lines=0
  if [[ -f "$tmp_out" ]]; then
    new_lines=$(wc -l < "$tmp_out" 2>/dev/null || echo 0)
  fi
  if [[ "$new_lines" -lt "$MIN_LINES" ]]; then
    echo "[WARN ] Incomplete local CSV ($new_lines lines): $tmp_out"
    echo "[FAIL clade] $tar_gz (incomplete output)" | tee -a "$FAIL_LOG"
    rm -rf "$dir_local" || true
    return 0
  fi

  mkdir -p "$parent"
  mv -f "$tmp_out" "$final_out"
  rm -rf "$dir_local" || true

  echo "[DONE ] CLADE  $final_out (lines=$new_lines)"
}
export -f process_one
export BASE REGIONS CLADE_RSCRIPT RSCRIPT_BIN MIN_LINES UNTAR_CHUNK COPY_LIMIT
export -f append_csv_no_dup_header _throttle _release_token

# ---------- discover and run ----------
mapfile -d '' TARS < <(find "$ROOT" -type f -name "hits_trees.tar.gz" -print0 | sort -z -V)
N_SRA="${#TARS[@]}"
if [[ "$N_SRA" -eq 0 ]]; then
  echo "No hits_trees.tar.gz under $ROOT. Nothing to do."
  exit 0
fi
if [[ "$N_SRA" -lt "$J" ]]; then J="$N_SRA"; fi
echo "Found $N_SRA SRAs; parallel jobs: $J"

printf '%s\0' "${TARS[@]}" \
  | xargs -0 -n1 -P "$J" bash -lc 'process_one "$@"' _

if [[ -s "$FAIL_LOG" ]]; then
  echo "Some SRAs failed in region $REGION:"
  cat "$FAIL_LOG"
else
  echo "All SRAs succeeded in region $REGION."
fi
