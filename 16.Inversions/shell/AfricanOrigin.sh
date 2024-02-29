## AFs for all Pops without missing data

mkdir /media/inter/mkapun/projects/DESTv2_data_paper/misc/Inversions/data

module load Tools/vcftools_0.1.13

## subset metadata and only retain samples with coverage >10
awk -F "," 'NR==1{print;next;} $11!="NA" && $7>10' /media/inter/mkapun/projects/ImPoolation/data/meta_cov.csv >/media/inter/mkapun/projects/DESTv2_data_paper/misc/Inversions/data/meta.csv

awk -F "," 'NR!=1 && $11!="NA" && $7>10 {print $1} ' /media/inter/mkapun/projects/ImPoolation/data/meta_cov.csv >/media/inter/mkapun/projects/DESTv2_data_paper/misc/Inversions/data/meta.txt

## convert to AF and only retain positions without missing data
gunzip -c /media/inter/mkapun/projects/ImPoolation/data/dest.all.PoolSNP.001.50.25Feb2023.norep.vcf.gz |
    awk '$0~/^\#/ || length($5)==1' |
    vcftools \
        --vcf - \
        --keep /media/inter/mkapun/projects/DESTv2_data_paper/misc/Inversions/data/meta.txt \
        --stdout \
        --recode |
    grep -v "\./\." |
    python /media/inter/mkapun/projects/ImPoolation/scripts/vcf2af.py \
        --input - |
    gzip >/media/inter/mkapun/projects/DESTv2_data_paper/misc/Inversions/data/DEST.af.gz

cp /media/inter/mkapun/projects/DESTv2_data_paper/misc/Inversions/data/DEST.af.gz tmp

## read breakpoints file
while IFS=$"," read -r inv chrom start end; do
    echo $inv
done </media/inter/mkapun/projects/DESTv2_data_paper/misc/Inversions/InvMeta/inversions_breakpoints_v5v6.txt

## split allele frequency dataset by genomic position of inversions
python /media/inter/mkapun/projects/DESTv2_data_paper/misc/Inversions/scripts/SplitByInv.py \
    --input /media/inter/mkapun/projects/DESTv2_data_paper/misc/Inversions/data/DEST.af.gz \
    --output /media/inter/mkapun/projects/DESTv2_data_paper/misc/Inversions/data/DEST_Inv \
    --param /media/inter/mkapun/projects/DESTv2_data_paper/misc/Inversions/InvMeta/inversions_breakpoints_v5v6.txt \
    --offset 200000

mkdir /media/inter/mkapun/projects/DESTv2_data_paper/misc/Inversions/results/AfricanOrigin

echo """
library(ggpubr)
library(tidyverse)
library("factoextra")
library("FactoMineR")
library(sp)

## read Meta data of reduced dataset (cov>10x)
Meta=read.csv("/media/inter/mkapun/projects/DESTv2_data_paper/misc/Inversions/data/meta.csv",
    header=T)

Meta$Sample <- Meta$sampleId
rownames(Meta)<-Meta$Sample

## read inversion frequency estimates for all samples
Inversion=read.table("/media/inter/mkapun/projects/DESTv2_data_paper/misc/Inversions/results/PoolSNP_nhm_inversion.af",
    header=T,
    na.string="na")   

## update colnames to match the naming scheme
colnames(Inversion)<-c("Sample","In2Lt","In2RNS","In3LP","In3RC","In3RK","In3RMo","In3RP")

## read AF data of genomic regions DISTANT to inversions
DATA.noinv=read.table(gzfile("/media/inter/mkapun/projects/DESTv2_data_paper/misc/Inversions/data/DEST_Inv_NoInv.af.gz"),
        header=T)  

## prepare AF matrix
AF.noinv=t(DATA.noinv[,3:ncol(DATA.noinv)])
rownames(AF.noinv)<-Meta$Sample

## make PCA
res.pca.noinv <- PCA(AF.noinv,  graph = FALSE)

### Plot results
Scree<-fviz_screeplot(res.pca.noinv, addlabels = TRUE, ylim = c(0, 50))
Scatter<-fviz_pca_ind(res.pca.noinv, label="none", habillage=as.factor(Meta$continent))
PLOT<-ggarrange(plotlist = list(Scree,Scatter),ncol=2,nrow=1,common.legend=TRUE, legend="bottom")
ggsave("/media/inter/mkapun/projects/DESTv2_data_paper/misc/Inversions/results/AfricanOrigin/PCA.pdf",
    PLOT,
    height=4,
    width=8)

### Get PC scores
ind.noinv <- get_pca_ind(res.pca.noinv)

## Now analyze all Inversions separately
plots <-  list()
plots2 <-  list()

for (i in c("In2Lt","In2RNS","In3LP","In3RP")) {
    print(i)

    ## read inv-specific AF file
    DATA=read.table(gzfile(paste0("/media/inter/mkapun/projects/DESTv2_data_paper/misc/Inversions/data/DEST_Inv_",i,".af.gz")),
        header=T)   

    ## prepare for PCA
    AF=t(DATA[,3:ncol(DATA)])
    rownames(AF)<-Meta$Sample

    ## do PCA and isolate PC scores
    res.pca <- PCA(AF,  graph = FALSE)
    ind <- get_pca_ind(res.pca)

    ## combine PC scores of first two axes from PCs based on AFs in inverted and non-inverted regions
    df<-data.frame(ind.noinv$coord[,1:2],ind$coord[,1:2],"Continent"=Meta$continent)
    colnames(df)<-c("PC1.noinv","PC2.noinv",paste0("PC1.",i),paste0("PC2.",i),"Continent")
    PC<-ggplot(df, aes_string(x=paste0("PC1.",i),y=paste0("PC2.",i),col="Continent"))+
        geom_point()+
        theme_bw()
    PC1<-ggplot(df, aes_string(x="PC1.noinv",y=paste0("PC1.",i),col="Continent"))+
        geom_point()+
        theme_bw()
    PC2<-ggplot(df, aes_string(x="PC2.noinv",y=paste0("PC2.",i),col="Continent"))+
        geom_point()+
        theme_bw()
    
    ## calulcate residuals of linear regressions between PC scores
    Res<-data.frame(
        "Continent"=Meta$continent,
        "Country"=Meta$country,
        "Sample"=Meta$Sample,
        "PC1.Res"=lm(ind.noinv$coord[,1]~ind$coord[,1])$residuals,
        "PC2.Res"=lm(ind.noinv$coord[,2]~ind$coord[,2])$residuals,
        "PC1"=ind$coord[,1],
        "PC2"=ind$coord[,2],
        "PC1.noinv"=ind.noinv$coord[,1],
        "PC2.noinv"=ind.noinv$coord[,2])
    
    ## plot
    Res.plot<-ggplot(Res, aes_string(x="PC1.Res",y="PC2.Res",col="Country"))+
        geom_point()+
        theme_bw()
    Plot1<-ggarrange(plotlist = list(PC,PC1,PC2,Res.plot),ncol=2,nrow=2,common.legend=TRUE, legend="bottom")
    FILE=paste0("/media/inter/mkapun/projects/DESTv2_data_paper/misc/Inversions/results/AfricanOrigin/PC_regressions_",i,".pdf")
    ggsave(FILE,
        Plot1,
        width=12,
        height=8)


    # Filter the points based on a condition
    filtered_points <- Res %>%
    filter(Continent=="Africa" & Country!="Morocco")  

    # Calculate the centroid of the filtered points
    centroid <- colMeans(filtered_points[, c("PC1.noinv", "PC2.noinv")])

    # Calculate the Euclidean distance between each point and the centroid
    Res$distances <- sqrt((Res$PC1.noinv - centroid[1])^2 + (Res$PC2.noinv - centroid[2])^2)

    Inv.t <- na.omit(Inversion %>%
        left_join(tibble::rownames_to_column(as.data.frame(ind$coord),var="Sample"), by ="Sample")) %>%
        left_join(Meta,by ="Sample") %>%
        left_join(Res,by="Sample") %>%
        select("Sample",i,"PC1.noinv","distances","country","continent")
    
    MINX=min(Inv.t[[i]])
    MINY=min(Inv.t$distances)+20
    fit1<-lm(subset(Inv.t, continent!="Africa" | country=="Marocco")$distance~subset(Inv.t, continent!="Africa" | country=="Marocco")[[i]])
    plots[[i]] <-ggplot(Inv.t,aes_string(x=i,y="distances",col="continent"))+
        geom_point()+
        ggtitle(i)+
        geom_smooth(data=subset(Inv.t, continent!="Africa" | country=="Marocco"),
            method = "lm", color = "black")+
        geom_label(aes(x=0,y = 200), hjust = 0,vjust=0, color="black",
             label = paste("Adj R2 = ",signif(summary(fit1)$adj.r.squared, 5),
                "\nIntercept =",signif(fit1$coef[[1]],5 ),
                " \nSlope =",signif(fit1$coef[[2]], 5),
                " \nP =",signif(summary(fit1)$coef[2,4], 5)),cex=2)+
        theme_bw()

    # Calculate the centroid of the filtered points
    centroid <- colMeans(filtered_points[, c("PC1", "PC2")])

    # Calculate the Euclidean distance between each point and the centroid
    Res$distances <- sqrt((Res$PC1 - centroid[1])^2 + (Res$PC2 - centroid[2])^2)

    Inv.t <- na.omit(Inversion %>%
        left_join(tibble::rownames_to_column(as.data.frame(ind$coord),var="Sample"), by ="Sample")) %>%
        left_join(Meta,by ="Sample") %>%
        left_join(Res,by="Sample") %>%
        select("Sample",i,"PC1","distances","country","continent")
    
    MINX=min(Inv.t[[i]])
    MINY=min(Inv.t$distances)+20
    fit1<-lm(subset(Inv.t, continent!="Africa" | country=="Marocco")$distance~subset(Inv.t, continent!="Africa" | country=="Marocco")[[i]])
    plots2[[i]] <-ggplot(Inv.t,aes_string(x=i,y="distances",col="continent"))+
        geom_point()+
        ggtitle(i)+
        geom_smooth(data=subset(Inv.t, continent!="Africa" | country=="Marocco"),
            method = "lm", color = "black")+
        geom_label(aes(x=0,y = 50), hjust=0,vjust=0, color="black",
             label = paste("Adj R2 = ",signif(summary(fit1)$adj.r.squared, 5),
                "\nIntercept =",signif(fit1$coef[[1]],5 ),
                " \nSlope =",signif(fit1$coef[[2]], 5),
                " \nP =",signif(summary(fit1)$coef[2,4], 5)),cex=2)+
        theme_bw()

}

Plot2<-ggarrange(plotlist = plots,ncol=2,nrow=2,common.legend=TRUE, legend="bottom")
FILE=paste0("/media/inter/mkapun/projects/DESTv2_data_paper/misc/Inversions/results/AfricanOrigin/final.pdf")
ggsave(FILE,
    Plot2,
    width=8,
    height=6)

Plot2.inv<-ggarrange(plotlist = plots2,ncol=2,nrow=2,common.legend=TRUE, legend="bottom")
FILE.inv=paste0("/media/inter/mkapun/projects/DESTv2_data_paper/misc/Inversions/results/AfricanOrigin/final_insideInv.pdf")
ggsave(FILE.inv,
    Plot2.inv,
    width=8,
    height=6)
