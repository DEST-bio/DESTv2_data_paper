#scp ~/ghcnd_records.RData aob2x@rivanna.hpc.virginia.edu:/project/berglandlab/alan/drosRTEC/mnt/pricey_1/dropPop/ghcnd_records.RData

### ijob -c 10 --mem=10G -p standard -A berglandlab
### module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### libraries
	library(data.table)
	library(foreach)
	library(doMC)
	registerDoMC(10)
  library(readxl)

### DrosRTEC samples
  samps <- fread("/scratch/aob2x/DESTv2/populationInfo/dest_v2.samps_13Jan2023.csv")

### GHCND data from Machado
  load("/project/berglandlab/alan/drosRTEC/mnt/pricey_1/dropPop/ghcnd_records.RData")
  gh <- as.data.table(ghcnd_records)

### conversion
  load("/project/berglandlab/alan/drosRTEC/mnt/pricey_1/dropPop/popSS.translate.Rdata")
  translate <- ss[,c("pop", "popOrig"), with=F]


### make sample conversion table
  gh.ag <- gh[,.N, list(Year, Locality=pop_name, pop=pop_name)]
  setkey(gh.ag, pop)
  setkey(translate, pop)

  gh.ag <- merge(gh.ag, translate)

### which samples were in the core20
  drosrtec <- as.data.table(read_excel("/scratch/aob2x/DESTv2/populationInfo/OriginalMetadata/elife-67577-supp1-v2.xlsx"))
  setkey(drosrtec, Locality, Year)
  gh.ag <- merge(gh.ag, drosrtec, by.x="popOrig", by.y="InternalName", all.y=T)

  tmp <- gh.ag[,c("Sample", "pop", "Core20"), with=F][Core20=="yes"]
  setkey(tmp, Sample)
  core <- tmp[!duplicated(tmp)]

### merge with new metadata
  core <- merge(core, samps, by.x="Sample", by.y="sampleId_orig")

### add new sample names back into GHCND
  gh <- merge(gh, core, by.x="pop_name", by.y="pop", allow.cartesian=T)

### load in NASA power data
  power.dt.ag <- foreach(fl.i=unique(gh$sampleId))%dopar%{
    # fl.i=unique(gh$sampleId)[1]
    message(fl.i)
    load(file=paste("~/DESTv2_data_paper/seasonal/nasapower_sampleId/", fl.i, ".Rdata", sep=""))
    power.dt.ag <- power.dt[,list(tmax.nasa=max(T2M), tmin.nasa=min(T2M)), list(sampleId, Year=YEAR, Month=MO, Day=DY)]
    power.dt.ag
  }
  power.dt.ag <- rbindlist(power.dt.ag)

### merge with gh data
  setkey(gh, sampleId, Year, Month, Day)
  setkey(power.dt.ag, sampleId, Year, Month, Day)

  m <- merge(gh, power.dt.ag)

  summary(lm(tmax~tmax.nasa*sampleId, m))

### save
  save(m, file="~/ghcnd_nasapower.Rdata")

### copy
  scp aob2x@rivanna.hpc.virginia.edu:~/ghcnd_nasapower.Rdata ~/.

  library(data.table)
  library(ggplot2)
  library(doMC)
  library(patchwork)
  library(viridis)
  registerDoMC(4)
  load("~/ghcnd_nasapower.Rdata")
  load("/Users/alanbergland/Documents/GitHub/dmel_seasonal_RTEC/predictability_model/weather_and_pbs.fix.Rdata")

### try to replicate the environmental models from
  m[,deltaDays:=as.numeric(make_date(Year, Month, Day) - (make_date(year) + jday - 1))]
  m[deltaDays>-14 & deltaDays<=0][,c("sampleId", "Year", "Month", "Day", "deltaDays")]

  m.use <- m[deltaDays>-14 & deltaDays<=0][,c("pop_name", "sampleId", "Year", "Month", "Day", "deltaDays", "sr_season", "tmax", "tmin", "tmax.nasa", "tmin.nasa", "locality", "exactDate")]
  m.use[,ly:=paste(locality, Year, sep="_")]

  #ggplot(data=m.use, aes(x=tmax, y=tmax.nasa)) + geom_point() + facet_wrap(~sampleId)

  m.use[,tmax.gh:=tmax/10]
  m.use[,tmin.gh:=tmin/10]

  o <- foreach(sp.th=seq(from=-50, to=400, by=5)/10, .combine="rbind")%dopar%{
    foreach(fall.th=seq(from=-50, to=400, by=5)/10, .combine="rbind")%do%{
      message(paste(sp.th, fall.th, sep=" / "))
      # sp.th <- 32; fall.th<-5
      m.ag <- m.use[exactDate==T,list(gh.tmax.spring=    mean(tmax.gh[sr_season=="spring"]>=sp.th, na.rm=T),      gh.tmin.fall=mean(tmin.gh[sr_season=="fall"]<=fall.th, na.rm=T),
                                      nasa.tmax.spring=mean(tmax.nasa[sr_season=="spring"]>=sp.th, na.rm=T), nasa.tmin.fall=mean(tmin.nasa[sr_season=="fall"]<=fall.th, na.rm=T)),
                                  list(ly, pop_name)]

      m.ag.ag <- merge(m.ag, pbs.small, by.x="pop_name", by.y="pop")
      gh.t1 <- summary(lm(beta~gh.tmax.spring + gh.tmin.fall, m.ag.ag))
      np.t1 <- summary(lm(beta~nasa.tmax.spring + nasa.tmin.fall, m.ag.ag))

      data.table(sp.th=sp.th, fall.th=fall.th,
                 r2=c(gh.t1$r.squared, np.t1$r.squared),
                 env=c("GHCND", "NASA"),
                 sp.cor=cor(m.ag.ag$gh.tmax.spring, m.ag.ag$nasa.tmax.spring, use="complete"),
                 fall.cor=cor(m.ag.ag$gh.tmin.fall, m.ag.ag$nasa.tmin.fall, use="complete"))

    }
  }


  ggplot(data=o, aes(x=sp.th, y=fall.th, fill=r2)) +
  geom_tile() +
  facet_grid(~env) + scale_fill_viridis_c()


  ggplot(data=o, aes(x=sp.th, y=sp.cor)) + geom_point()
  ggplot(data=o, aes(x=fall.th, y=fall.cor)) + geom_point()

  o[,list(sp=sp.th[which.max(r2)], fa=fall.th[which.max(r2)], r2=max(r2)), list(env)]

### mega plot
  temp.cor.tmax <- ggplot(data=m, aes(x=tmax.nasa, y=tmax/10)) + geom_point() +
  xlab("NASApower Tmax °C") + ylab("GHCND Tmax °C") +
  geom_abline(slope=1, intercept=0, color="red")

  temp.cor.tmin <- ggplot(data=m, aes(x=tmin.nasa, y=tmin/10)) + geom_point() +
  xlab("NASApower Tmin °C") + ylab("GHCND Tmin °C") +
  geom_abline(slope=1, intercept=0, color="red")

temp.cor.tmax + temp.cor.tmin


dest_v2[,date_string:=as.character(make_date(year) + jday - 1)]





### plays
m.ag <- m.use[exactDate==T,list(gh.tmax=max(tmax.gh, na.rm=T),
                                np.tmax=max(tmax.nasa, na.rm=T)),
                            list(ly, pop_name)]

m.ag.ag <- merge(m.ag, pbs.small, by.x="pop_name", by.y="pop")
gh.t1 <- summary(lm(beta~gh.tmax, m.ag.ag))
np.t1 <- summary(lm(beta~np.tmax, m.ag.ag))



  ggplot(data=m[sampleId=="US_Mas_Lan_1_2012-10-02"], aes(x=tmax, y=tmax.nasa)) + geom_point() + facet_grid(~sampleId)
