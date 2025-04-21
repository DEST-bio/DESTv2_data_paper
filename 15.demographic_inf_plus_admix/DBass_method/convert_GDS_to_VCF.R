# Install packages
#install.packages("BiocManager")
#install.packages("tidyverse")
#BiocManager::install("SeqArray")
#BiocManager::install("Rsamtools")

# This script trims and converts to VCF the GDS file in which the DEST2 data is
# stored. Running this script is the first step in analysis. Paths in the `Setup`
# section must be modified if used on a different machine.

# Setup
.libPaths("/home/dbass13/R/x86_64-conda-linux-gnu-library/4.2/")
library(SeqArray)
setwd("/data/rmccoy22/dbass13/DEST2/")
GDS_file <- "vcfs/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.gds"
metadata_file <- "metadata/dest_v2.samps_8Jun2023.csv"
VCF_file <- "vcfs/dest2_clustered_mainChroms.vcf.gz"

# Load the GDS File
GDS <- seqOpen(GDS_file)

# Use metadata to find and filter out samples that are not clustered
meta_df <- read.csv(metadata_file)
clustered_samps <- which(!is.na(meta_df$cluster2.0_k4))
seqSetFilter(GDS, sample.sel=clustered_samps)

# Filter out sex chromosomes and "chromosome" 4
seqSetFilterChrom(GDS, c("2L", "2R", "3L", "3R"))

# Output equivalent VCF file
seqGDS2VCF(GDS, VCF_file)
