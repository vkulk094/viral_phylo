# viral_phylo  
## Computational Workflow for Viral Phylogenetic and Diversity Analyses

This repository contains the full analysis pipeline developed for my master’s thesis project and the associated manuscript examining global patterns of viral biodiversity across latitudinal gradients. The repository is organized into clearly defined stages that progress from initial contig processing to phylogenetic analyses and final statistical visualizations.

All scripts are numbered to indicate execution order. Users should review the contents of each script before running, ensuring compatibility with their own computational environment, file paths, and software modules.

---

## Repository Structure

```
00_extract_cluster.py
00_submit_extract_clusters_array.sh

01_imgvr_cdhit_mafft_hmmer.R

02.1_N50_boxplots.py
02_contigs_assembly_megahit.sh

03_geNOMAD.sh

04_CheckV.sh

05_hmmsearch_after_mgm.sh

06_parse_hmmer.R
06_submit_parse_hmmer.sh

07_imgvr_uvig_phage_lookup_table...

08_Bridging_edge_IMGVR.R
08_submit_Bridging_edge_IMGVR.sh

09.1_IMGVR_boxplots_REBUILD.R
09_IMGVR_boxplots_v6.R

README.md
```

A brief description of each stage is provided below.

---

## Pipeline Overview

This pipeline performs the following major operations:

• Extracting clusters or sequence groups for downstream analysis  
• Running assembly and preprocessing steps  
• Annotating viral sequences  
• Performing quality assessments  
• Parsing HMM-based outputs  
• Integrating IMG/VR annotations  
• Calculating bridging edges or phylogenetic summaries  
• Producing final statistical and visual summary outputs  

A full schematic of the workflow will be added here once provided.

---

## Stage Descriptions

### 0. Cluster Extraction  
Files:  
`00_extract_cluster.py`  
`00_submit_extract_clusters_array.sh`  

This stage extracts cluster-level data required for downstream analysis. The included SLURM script defines array-based HPC execution parameters.

---

### 1. IMG/VR and Alignment Preparation  
File:  
`01_imgvr_cdhit_mafft_hmmer.R`  

Processes IMG/VR-associated sequences and prepares them for clustering, alignment, and similarity searches. All dependencies are handled within the script.

---

### 2. Assembly and Contig Length Summaries  
Files:  
`02_contigs_assembly_megahit.sh`  
`02.1_N50_boxplots.py`  

Performs metagenomic assembly and computes summary metrics such as N50. The Python script generates corresponding visual outputs for QC.

---

### 3. Viral Annotation  
File:  
`03_geNOMAD.sh`  

Executes annotation and classification steps in geNomad. Users should verify environment modules and adjust paths as required.

---

### 4. Quality Assessment  
File:  
`04_CheckV.sh`  

Runs viral quality assessment workflows on geNomad output.

---

### 5. Domain and Protein Search  
File:  
`05_hmmsearch_after_mgm.sh`  

Conducts HMM-based searches following gene prediction steps.

---

### 6. HMM Parsing  
Files:  
`06_parse_hmmer.R`  
`06_submit_parse_hmmer.sh`  

Parses HMM search results and formats outputs for subsequent analyses.

---

### 7. IMG/VR UVIG Phage Lookup Integration  
File:  
`07_imgvr_uvig_phage_lookup_tablemaker.R`  

Generates lookup tables linking IMG/VR UVIG records to phage-associated metadata.

---

### 8. Bridging Edge Computations  
Files:  
`08_Bridging_edge_IMGVR.R`  
`08_submit_Bridging_edge_IMGVR.sh`  

Computes bridging-edge metrics or related phylogenetic summaries.

---

### 9. Statistical Visualization (Final Figures)  
Files:  
`09.1_IMGVR_boxplots_REBUILD.R`  
`09_IMGVR_boxplots_v6.R`  

Produces final figure panels for the manuscript, including boxplots and summary statistical representations.
Note: The REBUILD file version allows for tweaking the diagrams without having to execute the computationally intensive steps that are prerequisite to figure construction.

---

## Dependencies

All dependencies invoked by these workflows are contained within the scripts.  
Users should verify compatibility with their own environment, adjusting:

• Executable paths  
• Software versions  
• Module or conda environment settings  
• Cluster-specific SLURM configurations  

---

## External Data Requirements

### IMG/VR v4.1 Database  
Download from:  
https://genome.jgi.doe.gov/portal/IMG_VR/IMG_VR.home.html  

Note: Registration with the JGI Genome Portal is required.  
The JGI Globus endpoint is recommended for transferring large datasets to HPC systems.

---

## Reproducibility Notes

• Execute scripts in numerical order.  
• SLURM scripts contain required HPC execution parameters.  
• Ensure that your environment provides all dependencies referenced within the scripts.  
• Final outputs include processed tables and publication-ready visualizations.

---

## Contact

For questions regarding workflow use or reproducibility:

**Vaibhav Kulkarni**  
University of Ottawa
Email: vkulk094@uottawa.ca
Website: https://vaibhavkulkarni.site
