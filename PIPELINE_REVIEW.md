# Pipeline summary and execution risk review

This document outlines the end-to-end viral phylogenetics workflow encoded in this repository and highlights potential execution pitfalls to double-check before running at scale.

## High-level pipeline flow
1. **Cluster extraction** (`00_extract_cluster.py`, `00_submit_extract_clusters_array.sh`): Extracts per-cluster protein FASTA files from an IMG/VR-wide index using cluster IDs provided at runtime. An accompanying SLURM array script parallelizes cluster retrieval.
2. **IMG/VR processing and alignment prep** (`01_imgvr_cdhit_mafft_hmmer.R`): Cleans IMG/VR sequences, clusters them, aligns with MAFFT, and prepares HMMER inputs.
3. **Assembly and contig metrics** (`02_contigs_assembly_megahit.sh`, `02.1_N50_boxplots.py`): Downloads SRA accessions, trims with Trimmomatic, assembles with MEGAHIT, and generates N50 QC plots.
4. **Annotation** (`03_geNOMAD.sh`): Runs geNomad annotation/classification on assembled contigs.
5. **Quality assessment** (`04_CheckV.sh`): Evaluates viral quality/completeness using CheckV.
6. **Domain/protein search** (`05_hmmsearch_after_mgm.sh`): Performs HMMER searches post-gene-prediction to annotate domains.
7. **HMM parsing** (`06_parse_hmmer.R`, `06_submit_parse_hmmer.sh`): Parses HMMER outputs and formats tables for downstream statistics.
8. **IMG/VR UVIG metadata integration** (`07_imgvr_uvig_phage_lookup_tablemaker.R`): Builds lookup tables linking IMG/VR UVIG records to phage metadata.
9. **Bridging-edge/phylogenetic metrics** (`08_Bridging_edge_IMGVR.R`, `08_submit_Bridging_edge_IMGVR.sh`): Computes graph-based bridging-edge metrics summarizing evolutionary relationships.
10. **Final statistics and figures** (`09.1_IMGVR_boxplots_REBUILD.R`, `09_IMGVR_boxplots_v6.R`): Reconstructs publication-ready boxplots and related summaries.

## Potential execution pitfalls to validate
- **Cluster extraction inputs and index consistency** (`00_extract_cluster.py`):
  - Assumes both `IMGVR_all_protein.faa` and its `all_viral.idx` index are present in the working directory and synchronized; rerunning after FASTA updates without rebuilding the index can return stale or missing entries.
  - Expects `cluster_seq_map.csv` to be comma-delimited with exactly two fields per line; malformed rows will raise `ValueError` and stop extraction mid-run because there is no error handling around parsing.
  - Command-line validation checks only the argument count; an empty or unknown cluster ID silently yields an empty FASTA.

- **Assembly SLURM orchestration** (`02_contigs_assembly_megahit.sh`):
  - Relies on a plain-text `sra_ids.txt` in the launch directory; absence immediately aborts the script. Duplicate IDs will overwrite checkpoints because job names and lockfiles reuse the raw ID.
  - The Trimmomatic adapter path is hardcoded (`/home/vaibhavk/adapters/TruSeq3-PE.fa`); running on another system without that file causes trimming to fail.
  - Resume jobs request three nodes but still set `-t 96` threads (less than the available cores), and the comment references two nodes—verify resource alignment with the target cluster scheduler.
  - Intermediate cleanup only removes `.fastq` and `.sra` files; partially generated MEGAHIT directories remain if a step fails before cleanup, and the lockfile is removed even on HMM/assembly errors, allowing concurrent reruns to clobber checkpoints.
  - Jobs are generated and immediately `sbatch`ed in a loop; failures in submission or downstream commands are not retried, and STDERR/STDOUT live inside per-job logs so global monitoring requires tailing each file.

- **Downstream R/ shell scripts** (stages 3–9):
  - Each script assumes required modules or packages are available (e.g., geNomad, CheckV, HMMER, R libraries) but does not perform dependency checks. Confirm module paths and versions before batch submission.
  - Many scripts expect prior directories and filenames exactly as produced by earlier steps; renaming outputs (e.g., moving MEGAHIT or CheckV results) will break subsequent relative-path assumptions.
  - SLURM submission scripts typically omit concurrency/array rate limits; on busy clusters, add `--array` throttling or `--requeue` flags to avoid oversubscription or silent preemption.

## Suggested preflight checks
- Verify presence and freshness of large reference files (IMG/VR FASTA, indexes, adapter sequences) on the intended filesystem.
- Test a single SRA ID end-to-end to validate module loads, scratch space quotas, and lockfile behavior before launching the full list.
- Centralize configuration (paths, thread counts, partition names) into a sourced file to minimize edits across numbered scripts.
- Add basic sanity checks (file existence, non-empty outputs) after each stage to prevent downstream scripts from operating on incomplete data.
