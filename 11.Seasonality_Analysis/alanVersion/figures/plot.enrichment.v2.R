
### plot
  library(ggplot2)
  library(data.table)

  system("scp aob2x@rivanna.hpc.virginia.edu://scratch/aob2x/DEST2_analysis/seasonality/glm_test_SEPT_29_2023/compiled/enrichment.NoCore20_seas_europe_yearPop_Ran.Rdata ~/.")
  load("~/enrichment.NoCore20_seas_europe_yearPop_Ran.Rdata")

### aggregate
  oo.ag <- oo[,list(N=sum(N), chr="genome"), list(perm, min_p, max_p, model_features, pops)]

  oo.ag.ag <- oo.ag[,list(N=N[perm==0],
                    N.perm.mean=mean(N[perm!=0]),
                    N.perm.lci=quantile(N[perm!=0], .005),
                    N.perm.uci=quantile(N[perm!=0], .995)),
                list(chr, min_p, max_p, pops, model_features)]

  save(oo.ag.ag, file="/Users/alanbergland/Documents/GitHub/DESTv2_data_paper/11.Seasonality_Analysis/alanVersion/figures/enrichment.NoCore20_seas_europe_yearPop_Ran.summary.Rdata")

  enrichment.plot <- ggplot(data=oo.ag.ag) +
  geom_ribbon(aes(x=log10(max_p), ymin=(N.perm.uci), ymax=N.perm.lci), alpha=.5, fill="grey") +
  geom_line(aes(x=log10(max_p), y=(N.perm.mean)), alpha=.5) +
  geom_line(aes(x=log10(max_p), y=(N)), color="red")

  ggsave(enrichment.plot, file="/Users/alanbergland/Documents/GitHub/DESTv2_data_paper/11.Seasonality_Analysis/alanVersion/figures.enrichment.NoCore20_seas_europe_yearPop_Ran.pdf")
