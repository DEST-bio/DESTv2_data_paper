# ijob -A berglandlab_standard -c4 -p dev --mem=20G
### module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R


.libPaths(c("/scratch/aob2x/Rlibs_4.3.1/")); .libPaths()

library(foreach)
library(data.table)
library(doMC)
registerDoMC(4)

fl <- list.files("/standard/vol186/bergland-lab/f3_admix/f3_full_array", full.name=T)

o <- foreach(fl.i=fl)%dopar%{
  message(fl.i); # fl.i <- fl[1]
  load(fl.i)

  as.data.table(f3.o.au)
}
o <- rbindlist(o)
setnames(o, "Z-score", "Z")
summary(o$f3)

summary(o$Z < -3)

o[,list(nSig=mean(Z < -3)), list(african_parent)][order(nSig)]

### local
  library(data.table)
  library(ggplot2)
  library(ggplot2)

  load("~/f3_megaoutput.Rdata")
  samps <- fread("/Users/alanbergland/Documents/GitHub/DESTv2/populationInfo/dest_v2.samps_8Jun2023.csv")
  setnames(samps, "cluster2.0_k4", "cluster")

  o <- merge(o, samps[,c("sampleId", "continent", "Recommendation", "lat", "locality")], by.x="focal", by.y="sampleId")
  o <- merge(o, samps[,c("sampleId", "continent", "Recommendation", "cluster", "country", "locality")], by.x="european_parent", by.y="sampleId")
  o <- merge(o, samps[,c("sampleId", "continent", "Recommendation", "cluster", "country")], by.x="african_parent", by.y="sampleId")

  setnames(o, c("Recommendation.x", "Recommendation.y", "Recommendation"), c("use_focal", "use_euro", "use_afr"))
  setnames(o, c("cluster.x", "cluster.y"), c("cluster_euro", "cluster_afr"))
  setnames(o, c("continent.x", "continent.y", "continent"), c("continent_focal", "continent_euro", "continent_afr"))
  setnames(o, c("country.x", "country.y"), c("country_euro", "country_afr"))
  setnames(o, c("locality.x", "locality.y"), c("locality", "locality_euro"))

  o <- o[use_focal=="Pass" & use_euro=="Pass" & use_afr=="Pass"]

  o[,cluster_euro:=as.character(cluster_euro)]
  o[,cluster_afr:=as.character(cluster_afr)]

  o[cluster_euro=="2", cluster_euro:="Eastern Europe"]
  o[cluster_euro=="3", cluster_euro:="Western Europe"]
  o[cluster_afr=="3", cluster_afr:="North Africa"]
  o[cluster_afr=="1", cluster_afr:="Sub-Saharan Africa"]

  o[country_euro=="UK", country_euro:="United Kingdom"]


#### first basic summary plot
  o[,p:=pnorm(Z, 0, 1)]
  o[,pa:=p.adjust(p, "bonferroni")]
  hist(o$p)
  table(o$pa<.05)

  o.ag <- o[,list(propSig=mean(pa<.05), Zmed=mean(Z), Zmed_sig=median(Z[pa<.05])), list(continent_focal, cluster_afr, country_euro, country_afr, cluster_euro)]
  o.ag[cluster_euro=="Western Europe", x:=1]
  o.ag[cluster_euro=="Eastern Europe", x:=.85]


    prop.sig <- ggplot(data=o.ag,
            aes(x=continent_focal,
                y=propSig,
                color=cluster_euro, group=interaction(cluster_euro, country_euro))) +
            geom_point() + geom_line() +
            geom_text(data=o.ag[cluster_afr=="Sub-Saharan Africa"][continent_focal=="North_America"],
              aes(x=x, y=propSig, label=interaction(country_euro)), size=2, hjust=1) +
            facet_grid(~cluster_afr+country_afr)

    Zmed.sig <- ggplot(data=o.ag[cluster_afr=="Sub-Saharan Africa"],
            aes(x=continent_focal,
                y=Zmed,
                color=cluster_euro, group=interaction(cluster_euro, country))) +
                geom_point() + geom_line() +
                geom_text(data=o.ag[cluster_afr=="Sub-Saharan Africa"][continent_focal=="North_America"],
                          aes(x=x, y=Zmed, label=interaction(country)), size=2, hjust=1)

    ZmedSig.sig <- ggplot(data=o.ag[cluster_afr=="Sub-Saharan Africa"],
            aes(x=continent_focal,
                y=Zmed_sig,
                color=cluster_euro, group=interaction(cluster_euro, country))) +
                geom_point() + geom_line() +
                geom_text(data=o.ag[cluster_afr=="Sub-Saharan Africa"][continent_focal=="North_America"],
                      aes(x=x, y=Zmed_sig, label=interaction(country)), size=2, hjust=1)



    prop.sig / Zmed.sig / ZmedSig.sig + plot_layout(guides="collect")

### Deeper dive into Germany
  hist(o[continent_focal=="North_America"][cluster_euro=="Eastern Europe"][country_euro=="Germany"]$Z)


  all_us <- ggplot(data=o[continent_focal=="North_America"][cluster_euro=="Eastern Europe"][country_euro=="Germany"][cluster_afr=="Sub-Saharan Africa"],
          aes(x=locality, y=Z, color=locality)) + geom_beeswarm() + coord_flip()

  cha_all_germany <- ggplot(data=o[continent_focal=="North_America"][country_euro=="Germany"][locality=="US_Vir_Cha"][cluster_afr=="Sub-Saharan Africa"][Z> -100],
          aes(x=focal, y=Z, color=locality_euro, group=interaction(focal, locality_euro))) + geom_boxplot() + coord_flip() + facet_grid(~cluster_afr)

all_us + cha_all_germany
ggplot(data=o[Z> -10000][continent_focal=="North_America"][cluster_afr=="Sub-Saharan Africa"], aes(x=lat, y=Z, group=lat)) + geom_boxplot()
