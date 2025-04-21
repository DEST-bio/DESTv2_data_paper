
require(poolfstat)
dmel.info=read.csv("/media/mathieu/DD12T_1/DIVERS_DATA_ANA/DEST_v2/dest_v2.samps_26April2023.csv",row.names=1)
dmel.full=readRDS("/media/mathieu/DD12T_1/DIVERS_DATA_ANA/DEST_v2/vcf/fullpooldataobj.Rds")

tmp.sim=which(dmel.full@poolnames=="SIM_SIM_w501_1_NA-MM-DD")
tmp.af=which(dmel.full@poolnames %in% c("ZA_Lim_Pha_1_2010-07-16","ZA_Mpu_Dul_1_2011-12-16","ZM_Sou_Sia_1_2010-07-16")) #purest according to Clustering (K5)
tmp.sel.eu=which(dmel.full@poolnames=="FR_Ill_Sai_1_2017-09-16")
tmp.sel.ena=which(dmel.info[dmel.full@poolnames,]$continent %in% c("North_America") & 
                    dmel.info[dmel.full@poolnames,]$lon> -100 &  
                    dmel.info[dmel.full@poolnames,]$Recommendation=="Pass")   
tmp.sel.wna=which(dmel.full@poolnames=="US_Cal_Sto_1_2013-09-03")

#coverage de simulans sont codes en 0 ou 1
dmel.sel=pooldata.subset(dmel.full,pool.index = c(tmp.sim,tmp.af,tmp.sel.eu,tmp.sel.wna,tmp.sel.ena),
                           min.maf = 0.01)
rm(dmel.full);gc()
tmp.sim=which(dmel.sel@poolnames=="SIM_SIM_w501_1_NA-MM-DD")
dmel.sel@readcoverage[,tmp.sim]=dmel.sel@readcoverage[,tmp.sim]*10
dmel.sel@refallele.readcount[,tmp.sim]=dmel.sel@refallele.readcount[,tmp.sim]*10
gc()


tmp.sel.ena=which(dmel.info[dmel.sel@poolnames,]$continent %in% c("North_America") & 
                  dmel.info[dmel.sel@poolnames,]$lon> -100 &  
                  dmel.info[dmel.sel@poolnames,]$Recommendation=="Pass") 
tmp.sel.ena=dmel.sel@poolnames[tmp.sel.ena]

for(i in 1:length(tmp.sel.ena)){
  tmp.sel=pooldata.subset(dmel.sel,pool.index = 
                          which(dmel.sel@poolnames%in%
                          c("FR_Ill_Sai_1_2017-09-16","SIM_SIM_w501_1_NA-MM-DD","ZA_Lim_Pha_1_2010-07-16",tmp.sel.ena[i],"US_Cal_Sto_1_2013-09-03")),
                          min.maf = 0.01,min.cov.per.pool = 5)
tmp.sel.f4=compute.fstats(tmp.sel,nsnp.per.bjack.block = 10000,computeF4 = T, return.F2.blockjackknife.samples = TRUE)
f4.ratio=compute.f4ratio(tmp.sel.f4,
                         num.quadruplet = paste0("FR_Ill_Sai_1_2017-09-16,SIM_SIM_w501_1_NA-MM-DD;ZA_Lim_Pha_1_2010-07-16,",tmp.sel.ena[i]),
                         den.quadruplet = "FR_Ill_Sai_1_2017-09-16,SIM_SIM_w501_1_NA-MM-DD;ZA_Lim_Pha_1_2010-07-16,US_Cal_Sto_1_2013-09-03")
if(i==1){res.f4.ratio=f4.ratio}else{res.f4.ratio=rbind(res.f4.ratio,f4.ratio)}
cat(tmp.sel.ena[i],"\n")

}

rownames(res.f4.ratio)=tmp.sel.ena
res.f4.ratio=cbind(res.f4.ratio,dmel.info[tmp.sel.ena,]$lat)
res.f4.ratio[,1:5]=1-res.f4.ratio[,1:5] #to estimate African ancestry
write.table(file="res.f4.ratio",res.f4.ratio,quote=F)

pdf("F4ratioEstimates.pdf",h=9,w=9)
plot(res.f4.ratio[,6],res.f4.ratio[,1],xlab="latitude",ylab="Prop. Of African Ancestry",pch=16,
     main="F4-ratio estimates of African ancestry in North-Eastern American samples")
text(res.f4.ratio[,6],res.f4.ratio[,1],tmp.sel.ena,cex=0.5)
graphics.off()
