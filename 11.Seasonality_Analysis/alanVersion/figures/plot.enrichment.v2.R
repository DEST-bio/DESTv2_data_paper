
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
  geom_ribbon(aes(x=-log10(max_p), ymin=(N.perm.uci), ymax=N.perm.lci), alpha=.5, fill="grey") +
  geom_line(aes(x=-log10(max_p), y=(N.perm.mean)), alpha=.5) +
  geom_line(aes(x=-log10(max_p), y=(N)), color="red") +
  xlab("-log10(p_lrt)")


  lab <- paste("Odds Ratio=", round(as.numeric(fisher.test(table(-log10(m$p_lrt)>3.5, m$C2_std_median>5))$estimate)), sep="")
  baypass_plot <- ggplot(m) +
  geom_hex(aes(x=-log10(p_lrt), y=C2_std_median)) +
  geom_vline(xintercept=3.3, color="red") +
  geom_hline(yintercept=5.5, color="red") + geom_text(aes(x=6, y=20, label="fuck you"))


  samps = fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_8Jun2023.csv")

  seasonal.sets <- samps[Recommendation=="Pass"][continent%in%c("Europe", "Asia")][set!="DrosRTEC"][,list(.N,
                                                                 season=c("spring", "fall"),
                                                                 sampleId=c(sampleId[which.min(jday)], sampleId[which.max(jday)]),
                                                                 jday=c(jday[which.min(jday)], jday[which.max(jday)]),
                                                                 diff=c(jday[which.min(jday)]-jday[which.max(jday)]),
                                                                 lat=lat[which.min(jday)]),
                                                              list(country, locality, year, loc_rep, cluster2.0_k5)][N>1][abs(diff)>30]

  collection_fig <- ggplot(data=seasonal.sets, aes(x=jday, y=lat, group=interaction(locality, year))) +
  geom_line() + ylab("Latitude") + xlab("Julian Day") + ggtitle("138 samples, 26 localities")


  seasonal_mega <- collection_fig + enrichment.plot + baypass_plot + plot_annotation(tag_level="A")

  ggsave(seasonal_mega, file="~/seasonal_composite.pdf")
  ggsave(enrichment.plot, file="/Users/alanbergland/Documents/GitHub/DESTv2_data_paper/11.Seasonality_Analysis/alanVersion/figures.enrichment.NoCore20_seas_europe_yearPop_Ran.pdf")
