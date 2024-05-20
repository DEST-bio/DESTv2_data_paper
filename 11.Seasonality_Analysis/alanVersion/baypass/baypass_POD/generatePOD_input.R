### arguments
  args <- commandArgs(trailingOnly=T)
  jobId <- as.numeric(args[1]); #jobId=51

  subpool <- jobId%%50; subpool
  set <- jobId - subpool; set

### libraries
  .libPaths(c("/scratch/aob2x/Rlibs_4.3.1/")); .libPaths()
  library(mvtnorm)

### source utility function
  source("/home/aob2x/baypass_afton/utils/baypass_utils.R")

### load allele freq distribution beta pramaters
  pbc1 <- read.table(paste0("/standard/vol186/bergland-lab/alan/dest_baypass/dest_subbaypass/destsubpool_", subpool, "_1_summary_beta_params.out"),h=T)$Mean
  pbc2 <- read.table(paste0("/standard/vol186/bergland-lab/alan/dest_baypass/dest_subbaypass/destsubpool_", subpool, "_2_summary_beta_params.out"),h=T)$Mean
  pbc3 <- read.table(paste0("/standard/vol186/bergland-lab/alan/dest_baypass/dest_subbaypass/destsubpool_", subpool, "_3_summary_beta_params.out"),h=T)$Mean
  pbc4 <- read.table(paste0("/standard/vol186/bergland-lab/alan/dest_baypass/dest_subbaypass/destsubpool_", subpool, "_4_summary_beta_params.out"),h=T)$Mean
  pbc5 <- read.table(paste0("/standard/vol186/bergland-lab/alan/dest_baypass/dest_subbaypass/destsubpool_", subpool, "_5_summary_beta_params.out"),h=T)$Mean

  pbc <- apply(do.call("rbind", list(pbc1, pbc2, pbc3, pbc4, pbc5)), 2, mean)

### load raw data
  data <- geno2YN(paste0("/standard/vol186/bergland-lab/alan/dest_baypass/dest_subpool/subpool_", subpool, ".genobaypass"))

### load omega mat
  omega1 <- read.table(paste0("/standard/vol186/bergland-lab/alan/dest_baypass/dest_subbaypass/destsubpool_", subpool, "_1_mat_omega.out"),h=F)
  omega2 <- read.table(paste0("/standard/vol186/bergland-lab/alan/dest_baypass/dest_subbaypass/destsubpool_", subpool, "_2_mat_omega.out"),h=F)
  omega3 <- read.table(paste0("/standard/vol186/bergland-lab/alan/dest_baypass/dest_subbaypass/destsubpool_", subpool, "_3_mat_omega.out"),h=F)
  omega4 <- read.table(paste0("/standard/vol186/bergland-lab/alan/dest_baypass/dest_subbaypass/destsubpool_", subpool, "_4_mat_omega.out"),h=F)
  omega5 <- read.table(paste0("/standard/vol186/bergland-lab/alan/dest_baypass/dest_subbaypass/destsubpool_", subpool, "_5_mat_omega.out"),h=F)

  omega <- as.matrix((omega1 + omega2 + omega3 + omega4 + omega5)/5)

### simulate
  setwd("/scratch/aob2x/dest2_baypass/pods_v2")
  tmp <- simulate.baypass(omega.mat = omega, nsnp = dim(data$NN)[1]*1, sample.size=data$NN, beta.pi=pbc, pi.maf=0, suffix=paste0("subpod_", jobId))
  write.table(omega, row.names=F, col.names=F, quote=F, sep="\t", file=paste0("omega.subpod_", jobId))
