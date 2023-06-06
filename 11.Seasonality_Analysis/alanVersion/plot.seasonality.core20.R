scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/DEST2_analysis/seasonality/compiled_output/enrichment.Core20_seas.Rdata ~/enrichment.Core20_seas.Rdata
scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/DEST2_analysis/seasonality/compiled_output/enrichment.Core20_seas.Rdata

system("scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/DEST2_analysis/seasonality/compiled_output/enrichment.NoCore20_seas.Rdata ~/.")



### libraries
  library(data.table)
  library(ggplot2)
  library(foreach)
  library(patchwork)

### data COre20
  load("~/enrichment.NoCore20_seas.Rdata")
  oo[,model_features:=factor(model_features, levels=c("LocBinomial", "LocQB", "PhyloQB", "Loc_PhyloQB", "LocRan", "Phylo_LocRan"))]

  oo.ag <- oo[,list(N=sum(N), chr="genome"), list(perm, max_p, model_features, pops)]

  oo <- rbind(oo, oo.ag, fill=T)

####
  oo.ag <- oo[,list(en=median(log2(N[perm==0]/N[perm!=0]))), list(pops, max_p, model_features)]

  p1 <-
  ggplot(data=oo[perm!=0][chr=="genome"]) +
  geom_line(aes(x=log10(max_p), y=(N), group=perm), alpha=.5) +
  geom_line(data=oo[perm==0][chr=="genome"], aes(x=log10(max_p), y=(N), group=perm), color="red") +
  facet_grid(model_features~chr, scales="free_y")

  p2 <-
  ggplot(data=oo[perm!=0][chr!="genome"]) +
  geom_line(aes(x=log10(max_p), y=(N), group=perm), alpha=.5) +
  geom_line(data=oo[perm==0][chr!="genome"], aes(x=log10(max_p), y=(N), group=perm), color="red") +
  facet_grid(model_features~chr, scales="free_y")

layout <-
"ABBBB"

  mega <- p1 + p2 + plot_layout(design=layout)
  ggsave(mega, file="~/NoCore20.png", h=10, w=10)


  ggplot(data=oo.ag) +
  geom_line(aes(x=log10(max_p), y=en), alpha=.5) +
  facet_wrap(~model_features, scales="free")
