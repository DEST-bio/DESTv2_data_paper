# ijob -A berglandlab_standard -c5 -p standard --mem=4G
### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### libraries
  library(data.table)
  library(foreach)

### compile data
  fl <- list.files("/scratch/yey2sn/DEST2_analysis/seasonality/GLM_omnibus_JULY_19_2023/compiled/enrichment", full.names=T)

  oo <- foreach(fl.i=fl, .combine="rbind")%do%{
    # fl.i <- fl[1]
    message(fl.i)
    load(fl.i)
    oo
  }

  save(oo, file="/scratch/aob2x/DEST2_analysis/seasonality/GLM_omnibus_JUNE_5_2023/compiled/mega_enrichment.JUNE_13_2023.Rdata")

### plot
  library(ggplot2)
  library(data.table)

  system("scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/DEST2_analysis/seasonality/GLM_omnibus_JUNE_13_2023/compiled/mega_enrichment.JUNE_5_2023.Rdata ~/.")
  load("~/mega_enrichment.JUNE_13_2023.Rdata")

### collapse

  oo <- oo[,list(N=sum(N)), list(min_p, max_p, perm, pops, model_features, chr)]

  oo[,model_features:=factor(model_features,
                      levels=c("yearPop_Binomial", "Phylo_yearPop_Binomial",      "yearPop_Ran",      "Phylo_yearPop_Ran",
                                "Loc_Binomial_PCA",  "Loc_Ran_PCA"))]

  oo.ag <- oo[,list(N=sum(N), chr="genome"), list(perm, max_p, model_features, pops)]

  oo <- rbind(oo, oo.ag, fill=T)

####
  foreach(pops.i=unique(oo$pops))%do%{
    message(pops.i)

    p1 <-
    ggplot(data=oo[perm!=0][chr=="genome"][pops==pops.i]) +
    geom_line(aes(x=log10(max_p),
    y=(N),
    group=perm),
     alpha=.5) +
    geom_line(data=oo[perm==0][chr=="genome"][pops==pops.i],
     aes(x=log10(max_p),
     y=(N),
     group=perm),
     color="red") +
     facet_grid(#model_features
     ~chr, scales="free_y")
ggsave(p1,file = "p1.pdf")

    p2 <-
    ggplot(data=oo[perm!=0][chr!="genome"][pops==pops.i]) +
    geom_line(aes(x=log10(max_p), y=(N), group=perm), alpha=.5) +
    geom_line(data=oo[perm==0][chr!="genome"][pops==pops.i], aes(x=log10(max_p), y=(N), group=perm), color="red") +
    facet_grid(#model_features
    ~chr, scales="free_y")
ggsave(p2,file = "p2.pdf")

    layout <-
    "ABBBB"

    mega <- p1 + p2 + plot_layout(design=layout) +  plot_annotation(title = pops.i)
<<<<<<< HEAD
    ggsave(mega, file=paste("~/modelEnrichment.JUNE_13_2023.", pops.i, ".png", sep=""), h=8, w=8)
=======

    ggsave(mega, file=paste("./modelEnrichment.new.", pops.i, ".png", sep=""), h=8, w=8)
>>>>>>> 7ee24b1bdddc9e0b332daa88c532dda6e6d27502

  }
