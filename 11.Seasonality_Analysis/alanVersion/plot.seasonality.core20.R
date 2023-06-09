scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/DEST2_analysis/seasonality/compiled_output/enrichment.Core20_seas.Rdata ~/enrichment.Core20_seas.Rdata
scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/DEST2_analysis/seasonality/compiled_output/enrichment.Core20_seas.Rdata

system("scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/DEST2_analysis/seasonality/compiled_output/enrichment.NoCore20_seas.Rdata ~/.")
system("scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/DEST2_analysis/seasonality/compiled_output/enrichment.NoCore20_NoProblems_seas.Rdata ~/.")
system("scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/DEST2_analysis/seasonality/compiled_output/enrichment.all.Rdata ~/.")

### libraries
  library(data.table)
  library(ggplot2)
  library(foreach)
  library(patchwork)

### data COre20
  load("~/enrichment.all.Rdata")
  oo[,model_features:=factor(model_features, levels=c("LocBinomial", "LocQB", "PhyloQB", "Loc_PhyloQB", "LocRan", "Phylo_LocRan"))]

  oo.ag <- oo[,list(N=sum(N), chr="genome"), list(perm, max_p, model_features, pops)]

  oo <- rbind(oo, oo.ag, fill=T)

####
  foreach(pops.i=unique(oo$pops))%do%{
    message(pops.i)

    p1 <-
    ggplot(data=oo[perm!=0][chr=="genome"][pops==pops.i]) +
    geom_line(aes(x=log10(max_p), y=(N), group=perm), alpha=.5) +
    geom_line(data=oo[perm==0][chr=="genome"][pops==pops.i], aes(x=log10(max_p), y=(N), group=perm), color="red") +
    facet_grid(model_features~chr, scales="free_y")

    p2 <-
    ggplot(data=oo[perm!=0][chr!="genome"][pops==pops.i]) +
    geom_line(aes(x=log10(max_p), y=(N), group=perm), alpha=.5) +
    geom_line(data=oo[perm==0][chr!="genome"][pops==pops.i], aes(x=log10(max_p), y=(N), group=perm), color="red") +
    facet_grid(inv+model_features~chr, scales="free_y")

    layout <-
    "ABBBB"

    mega <- p1 + p2 + plot_layout(design=layout) +  plot_annotation(title = pops.i)
    ggsave(mega, file=paste("~/modelEnrichment.", pops.i, ".png", sep=""), h=5, w=10)

  }



  ggplot(data=oo.ag) +
  geom_line(aes(x=log10(max_p), y=en), alpha=.5) +
  facet_wrap(~model_features, scales="free")
