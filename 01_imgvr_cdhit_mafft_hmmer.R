# dependencies: 
# - via pip install: seqmagick (permanent install for user, python)
# - via brew install:  maftt, hmmer3

# module load StdEnv/2020 r/4.2.2 python/3.10.2 mafft/7.471 hmmer/3.3.2 

# install.packages(c("foreach", "doParallel", "seqinr"))

library(seqinr)                 # Load relevant libraries
library(doParallel)
library(foreach)

# Detect Available Cores
available_cores <- parallel::detectCores(logical = TRUE)
registerDoParallel(cores = available_cores)
message(paste("Detected", available_cores, "cores. Using all available for parallel processing."))

# Parse CD-HIT Output
infile <- readLines("all_viral_cd_hit.clstr")  # Read the CD-HIT output file
nLines <- length(infile)

# Initialize Variables
cluster_lst <- list()  # List to store sequences for each cluster
cluster_name <- c()    # Current cluster name
cluster_names <- c()   # Names of all clusters

# Parse clusters
for (i in 1:nLines) {
    if (grepl("Cluster", infile[i])) {
        cluster_name <- sub(">", "", infile[i])
        cluster_name <- sub(" ", "_", cluster_name)
        cluster_names <- c(cluster_names, cluster_name)
        print(paste0("Now parsing cluster: ", cluster_name))
    } else if (length(strsplit(infile[i], " ")[[1]]) > 1) {
        tmp1 <- strsplit(infile[i], " ")
        tmp2 <- sub(">", "", unlist(tmp1)[2])  # Extract sequence ID
        tmp3 <- sub("\\.\\.\\.", "", tmp2)     # Remove trailing "..."
        cluster_lst[[cluster_name]] <- c(cluster_lst[[cluster_name]], tmp3)
    } else {
        message(paste("Unexpected line format:", infile[i]))
    }
}

# Debug Cluster Parsing
print(head(cluster_lst))  # Check the first few clusters

# Read and Clean Reference FASTA
refseq_faa <- read.fasta("all_viral.faa")  # Load the reference FASTA file

# Extract only base IDs (remove descriptors like protein names)
names(refseq_faa) <- sub("^>(\\S+).*", "\\1", names(refseq_faa))  # Keep only the first word in the header

# Debug: Check the cleaned sequence names
print(head(names(refseq_faa)))

# Create Output Directories
dir.create("clusters", showWarnings = FALSE)
dir.create("clusters/faa", showWarnings = FALSE)
dir.create("clusters/aln", showWarnings = FALSE)

# Define Minimum Sequence Threshold
min_nb_seq <- 10  # Minimum number of sequences for a cluster to be processed

# Process Clusters and Write .faa Files
cluster_names2aln <- c()  # List of cluster names to align

for (i in cluster_names) {
    print(paste0("Now processing cluster: ", i, " out of ", length(cluster_names)))

    acc_lst <- cluster_lst[[i]]
    if (length(acc_lst) < min_nb_seq) {
        message(paste("Skipping cluster", i, "due to insufficient sequences"))
        next
    }

    # Match cleaned IDs in the reference FASTA
    pos <- which(names(refseq_faa) %in% acc_lst)
    if (length(pos) == 0) {
        message(paste("No matching sequences found in refseq_faa for cluster:", i))
        message(paste("Cluster IDs:", paste(acc_lst, collapse = ", ")))
        next
    }

    tmp_seq <- refseq_faa[pos]
    write.fasta(tmp_seq, acc_lst, paste0("clusters/faa/", i, ".faa"))
    cluster_names2aln <- c(cluster_names2aln, i)
    print(paste("Cluster", i, "written to file."))
}

# Align Sequences and Build HMMs
dir.create("clusters/sto", showWarnings = FALSE)
dir.create("clusters/hmm", showWarnings = FALSE)

aln_created <- foreach(i = 1:length(cluster_names2aln), .combine = rbind) %dopar% {
    cluster_name <- cluster_names2aln[i]
    faa_file <- paste0("clusters/faa/", cluster_name, ".faa")
    aln_file <- paste0("clusters/aln/", cluster_name, ".fasta")
    sto_file <- paste0("clusters/sto/", cluster_name, ".sto")
    hmm_file <- paste0("clusters/hmm/", cluster_name, ".hmm")

    # Debug: Check if .faa file exists and is non-empty
    if (!file.exists(faa_file) || file.size(faa_file) == 0) {
        message(paste("Skipping empty or missing FAA file:", faa_file))
        return(NA)
    }

    # Run MAFFT
    mafft_status <- system(paste0("mafft --auto ", faa_file, " > ", aln_file))
    if (mafft_status != 0 || !file.exists(aln_file) || file.size(aln_file) == 0) {
        message(paste("MAFFT failed for file:", faa_file))
        return(NA)
    }

    # Convert to Stockholm format
    seqmagick_status <- system(paste0("seqmagick convert ", aln_file, " ", sto_file))
    if (seqmagick_status != 0 || !file.exists(sto_file) || file.size(sto_file) == 0) {
        message(paste("Seqmagick failed for file:", aln_file))
        return(NA)
    }

    # Build HMM
    hmmbuild_status <- system(paste0("hmmbuild --amino --cpu 4 ", hmm_file, " ", sto_file))
    if (hmmbuild_status != 0 || !file.exists(hmm_file) || file.size(hmm_file) == 0) {
        message(paste("HMMER failed for file:", sto_file))
        return(NA)
    }

    return(cluster_name)
}

# Combine HMM Files
system("cat clusters/hmm/*.hmm > clusters/vir_models.hmm")

# Exit R Session
q(save = "no")
