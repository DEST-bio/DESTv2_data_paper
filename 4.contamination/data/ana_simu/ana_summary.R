n.ind.simu=seq(0,80,2)
n.ind.mela=80-n.ind.simu

nkmin=c(1,5) ; conf=c(0.90,0.95)
cdt=expand.grid(nkmin,conf)

res=array(0,dim=c(length(n.ind.simu),length(nkmin)*length(conf),6),
            dimnames = list(n.ind.simu,paste0("Nmin",cdt[,1]," conf",cdt[,2]),
          c("NseqSimu","NseqMela","NseqAnan","NseqWolb","NseqOthers","NseqAssigned")))

for(i in 1:length(n.ind.simu)){
  for(j in 1:nrow(cdt)){
    tmp.f=paste0("res_files/sim.",n.ind.simu[i],"_nkmin",cdt[j,1],"_conf",cdt[j,2],".summary")
    tmp.d=read.table(tmp.f,skip=1,row.names=1)
    tmp.v=c(tmp.d["Dsimu",1],tmp.d["Dmela",1],tmp.d["Danan",1],tmp.d["Wolb",1])
    res[i,j,]=c(tmp.v,sum(tmp.d[,1])-sum(tmp.v),sum(tmp.d[,1]))
  }
}

pdf("res_simu.pdf",h=18,w=18)
tmp.layout=c(1,1,2,3)
#tmp.layout=rbind(tmp.layout,tmp.layout,3:6)
tmp.layout=cbind(tmp.layout,c(4,5,6,6))
layout(tmp.layout)

#Relative percentage of Dsimu
tmp.d=100*res[,,1]/(res[,,1]+res[,,2])
rmse=sqrt(colSums((tmp.d/100-n.ind.simu/80)**2)/nrow(tmp.d))
# Nmin1 conf0.9  Nmin5 conf0.9 Nmin1 conf0.95 Nmin5 conf0.95
# 0.05222908     0.05531143     0.04926308     0.05207857
matplot(n.ind.simu,tmp.d,xlab="Number of D. simulans indiv. in simulation",
        ylab="%",pch=16,type="b",lty=2,main="A) Relative Proportion of simulans vs. melanogaster seq")
points(n.ind.simu,100*n.ind.simu/80,type="l",col="grey",lty=1)
legend("topleft",paste0("Nk>",cdt[,1]," C>",cdt[,2]," (RMSE=",round(rmse,digits=3),")"),col=1:4,pch=16,lty=2,ncol=2)
legend("bottomright","Exp. Prop",col="grey",lty=1,bty="n")

#percentage of seq assigned
tmp.d=res[,,6]
matplot(n.ind.simu,100*tmp.d/80000,xlab="Number of D. simulans indiv. in simulation",
        ylab="%Assigned Seq",pch=16,type="b",lty=2,main="B) Percentage of Assigned Reads")
legend("topright",paste0("Nk>",cdt[,1]," C>",cdt[,2]),col=1:4,pch=16,lty=2,ncol=2)

#percentage of sequences assigned to other species than melano or simulans
tmp.d=100*(res[,,6]-res[,,1]-res[,,2])/res[,,6]
matplot(n.ind.simu,tmp.d,xlab="Number of D. simulans indiv. in simulation",
        ylab="%",pch=16,type="b",lty=2,main="C) Percentage of sequences neither assigned to Dsimu nor Dmela")
legend("topleft",paste0("Nk>",cdt[,1]," C>",cdt[,2]),col=1:4,pch=16,lty=2,ncol=2)

#percentage of Wolbachia
tmp.d=100*res[,,4]/res[,,6]
matplot(n.ind.simu,tmp.d,xlab="Number of D. simulans indiv. in simulation",
        ylab="%",pch=16,type="b",lty=2,main="D) Percentage of sequences assigned to Wolbachia")
legend("topleft",paste0("Nk>",cdt[,1]," C>",cdt[,2]),col=1:4,pch=16,lty=2,ncol=2)

#percentage of Dananassae
tmp.d=100*res[,,3]/res[,,6]
matplot(n.ind.simu,tmp.d,xlab="Number of D. simulans indiv. in simulation",
        ylab="%",pch=16,type="b",lty=2,main="E) Percentage of sequences assigned to D. ananassae")
legend("top",paste0("Nk>",cdt[,1]," C>",cdt[,2]),col=1:4,pch=16,lty=2,ncol=2)

#Relative percentage of Dsimu (with emp. correction for uneven Wolbachia contaminant using estimates on data with n=0 and n=80)
ps=as.matrix(n.ind.simu/80)
tmp.d=100*(res[,,1]+ps%*%res[41,,4])/(res[,,1]+res[,,2]+ps%*%res[41,,4]+(1-ps)%*%res[1,,4])
rmse=sqrt(colSums((tmp.d/100-n.ind.simu/80)**2)/nrow(tmp.d))
# Nmin1 conf0.9  Nmin5 conf0.9 Nmin1 conf0.95 Nmin5 conf0.95
# 0.01669223     0.01538168     0.01419915     0.01177608
matplot(n.ind.simu,tmp.d,xlab="Number of D. simulans indiv. in simulation",
        ylab="%",pch=16,type="b",lty=2,main="F) Relative Proportion of simulans vs. melanogaster seq (with emp. correction for unequal Wolbachia contamination)")
points(n.ind.simu,100*n.ind.simu/80,type="l",col="grey",lty=1)
legend("topleft",paste0("Nk>",cdt[,1]," C>",cdt[,2]," (RMSE=",round(rmse,digits=3),")"),col=1:4,pch=16,lty=2,ncol=2)
legend("bottomright","Exp. Prop",col="grey",lty=1,bty="n")

graphics.off()
