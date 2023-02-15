### load in the datasets 
### 
library(vroom)
library(tidyverse)
###
seqtk <- "/home/yey2sn/software/seqtk/seqtk sample"
###
sim.samps <- vroom("/scratch/yey2sn/DEST2_analysis/kmer.analysis/sim_fasta/Dsim.80.txt")
mel.samps <- vroom("/scratch/yey2sn/DEST2_analysis/kmer.analysis/mel_fasta/Dmel.80.txt")

####
####
args = commandArgs(trailingOnly=TRUE)
k=as.numeric(args[1])

##foreach(k = 0:40, .combine = "rbind")%do%{
  
  sim.share = 80*(k/(40))
  mel.share = 80-(80*(k/(40)))
  
  ##
  message(paste(mel.share, sim.share, sim.share+mel.share , sep = "/"))
  ##
  
  ##
  if(sim.share != 0){
    sim.samps[1:sim.share,"Run"] -> chosen.sims
  }
  
  ## 
  if(mel.share != 0){
    mel.samps[1:mel.share,"Run"] -> chosen.mels
  }
  ##
  
  out.name = paste("sim",sim.share, "mel",mel.share, "tot",sim.share+mel.share, sep = "." )
  system(paste("rm", paste(out.name, ".1000r.fq", sep =""), sep = " "))
  
  if(sim.share != 0){
  message("pass sim")
  for(i in 1:dim(chosen.sims)[1]){
    message(paste("now processing", "sim" , chosen.sims$Run[i]))
    system(paste(seqtk, paste("./sim_fasta/",chosen.sims$Run[i],".fastq.gz", sep = "") ,
                 "1000", ">>" , 
                 paste(out.name, ".1000r.fq", sep =""), 
                 sep = " " ))   
  }
  } ## dim check!
  
  if(mel.share != 0){
    message("pass mel")
  for(j in 1:dim(chosen.mels)[1]){
    message(paste("now processing", "mel" , chosen.mels$Run[j]))
    system(paste(seqtk, paste("./mel_fasta/",chosen.mels$Run[j],".fastq.gz", sep = "") ,
                 "1000", ">>" , 
                 paste(out.name, ".1000r.fq", sep =""), 
                 sep = " " ))   
  }
  } ## mel check
  
###}


