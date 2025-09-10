#!/bin/bash
#SBATCH --job-name=extract_clust_batch
#SBATCH --output=/home/vaibhavk/scratch/IMG_VR/logs/extract_clust_batch_%A_%a.out
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=02:59:00    # ~193 clusters at ~0.75min/cluster = ~2.4h + buffer
#SBATCH --array=1-9957     # 1,921,454 total / 193 per job ≃ 9957 tasks, up to 500 concurrent


module load StdEnv/2020 python/3.11.2
cd /home/vaibhavk/scratch/IMG_VR

BATCH_SIZE=193
TOTAL=$(wc -l < cluster_ids.txt)

# Compute the slice for this array task
START=$(( (SLURM_ARRAY_TASK_ID-1)*BATCH_SIZE + 1 ))
END=$(( SLURM_ARRAY_TASK_ID*BATCH_SIZE ))
[ $END -gt $TOTAL ] && END=$TOTAL

# Run the extractor on each cluster in this batch
sed -n "${START},${END}p" cluster_ids.txt | while read CLUSTER_ID; do
    ./extract_cluster.py "$CLUSTER_ID"
done