scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/DEST2_analysis/seasonality/compiled_output/enrichment.Core20_seas.Rdata ~/enrichment.Core20_seas.Rdata
scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/DEST2_analysis/seasonality/compiled_output/enrichment.Core20_seas.Rdata

scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/DEST2_analysis/seasonality/compiled_output/enrichment.NoCore20_seas.Rdata ~/.



### libraries
  library(data.table)
  library(ggplot2)
  library(foreach)

### data
  load("~/enrichment.NoCore20_seas.Rdata")

####
  oo.ag <- oo[,list(en=median(log2(N[perm==0]/N[perm!=0]))), list(pops, max_p, model_features)]

  ggplot(data=oo[perm!=0]) +
  geom_line(aes(x=log10(max_p), y=(N), group=perm), alpha=.5) +
  geom_line(data=oo[perm==0], aes(x=log10(max_p), y=(N), group=perm), color="red") +
  facet_wrap(~model_features, scales="free")


  ggplot(data=oo.ag) +
  geom_line(aes(x=log10(max_p), y=en), alpha=.5) +
  facet_wrap(~model_features, scales="free")
