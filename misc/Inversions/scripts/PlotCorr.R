SNAPE=read.table("/media/DEST2_NHM/results/SNAPE_nhm_inversion.af",header=T)
PoolSNP=read.table("/media/DEST2_NHM/results/PoolSNP_nhm_inversion.af",header=T)
PoolSNP2 <- PoolSNP[PoolSNP$Sample %in% SNAPE$Sample,] 
header=colnames(SNAPE)[2:ncol(SNAPE)]
header=header[header!="In.3R.K"]
pdf("/media/DEST2_NHM/results/SNAPEvsPoolSNP.pdf",width=12,height=8)
par(mfrow=c(2,3))
for (i in header){
  print(i)
  reg=lm(as.numeric(paste(SNAPE[[i]]))~as.numeric(paste(PoolSNP2[[i]])))
  plot(as.numeric(paste(PoolSNP2[[i]])),as.numeric(paste(SNAPE[[i]])),xlab="SNAPE",ylab="PoolSNP",pch=16,col=rgb(0,0,0,0.2),main=i)
  abline(reg,col="red",lty=2,lwd=2)
  #legend("topleft",legend=substitute(expression(italic(R)^2 == myvalue),myvalue =round(summary(reg)$'r.squared',2)),box.col="NA")
}
dev.off()
