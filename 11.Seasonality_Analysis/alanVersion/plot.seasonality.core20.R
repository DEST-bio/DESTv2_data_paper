scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/DEST2_analysis/seasonality/compiled_output/enrichment.Core20_seas.Rdata ~/enrichment.Core20_seas.Rdata

### libraries
  library(data.table)
  library(ggplot2)
  library(foreach)

### data
  load("~/enrichment.Core20_seas.Rdata")

####


  ggplot(data=oo[perm!=0]) +
  geom_line(aes(x=log10(thr), y=(nSig), group=perm), alpha=.5) +
  geom_line(data=oo[perm==0], aes(x=log10(thr), y=(nSig), group=perm), color="red") +
  facet_wrap(~model_features)
