## get SYNC data
scp mkapun@10.95.0.14:/media/DEST2_NHM/output/* /media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/data
scp mkapun@10.95.0.14:/media/DEST2_NHM/data/SNAPE.header.txt /media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/data/

mkdir /media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/data

## PoolSNP
gunzip -c /media/inter/ssteindl/DEST/DEST2_NHM/collapsed/PoolSNP/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.vcf.gz |
    parallel \
        --jobs 200 \
        --pipe \
        -k \
        --cat python3 /media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/scripts/VCF2sync.py \
        --input {} |
    gzip >/media/inter/ssteindl/DEST/DEST2_NHM/collapsed/PoolSNP/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.sync.gz

## recode missing data confoming with Grenedalf .:.:.:.:.:. > 0:0:0:0:0:0
pigz -dc /media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/data/dest.PoolSeq.SNAPE.NA.NA.25Feb2023.norep.gz |
    sed 's/.:.:.:.:.:./0:0:0:0:0:0/g' |
    pigz >/media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/data/DEST2_snape.sync.gz

pigz -dc /media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/data/dest.all.PoolSNP.001.50.25Feb2023.norep.gz |
    sed 's/.:.:.:.:.:./0:0:0:0:0:0/g' |
    pigz >/media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/data/DEST2_poolsnp.sync.gz

pigz -dc /media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/data/dest.all.PoolSNP.001.50.25Feb2023.norep.gz |
    sed 's/.:.:.:.:.:./0:0:0:0:0:0/g' |
    pigz >/media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/data/DEST2_poolsnp.sync.gz

## Define populations with poolsize < 15 which causes problems in Grenedalf when the populations have also low coverage
echo '''CM_Nor_Oku_1_2004-04-15
ET_Amh_Deb_1_2011-12-16
ET_Oro_Dod_1_2008-12-16
ET_Sou_Bon_1_2011-12-16
GA_Hau_Fra_1_2002-03-16
GN_Dal_Don_1_2005-06-15
KE_Cen_Nya_1_2009-01-16
SIM_SIM_w501_1_NA-MM-DD
NG_Bor_Mai_1_2004-09-16
UG_Kis_Kis_1_2012-01-16
US_Vir_Cha_1_2018-06-28
ZA_Eas_Bar_1_2011-12-16
ZW_Mat_Vic_1_2001-07-16''' >/media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/data/exclude.txt

## Define Chromsomal arms to consider
Chrom=("4" "2L" "2R" "3L" "3R" "X")
Lengths=(1348131 23513712 25286936 28110227 32079331 23542271)
threads=150

## exclude populations
for i in poolsnp snape; do

    echo """
    #!/bin/sh

    ## name of Job
    #PBS -N mapping

    ## Redirect output stream to this file.
    #PBS -o /media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/results/${i}_filterSync_log.txt

    ## Stream Standard Output AND Standard Error to outputfile (see above)
    #PBS -j oe

    ## Select a maximum of ${threads} cores and 800gb of RAM
    #PBS -l select=1:ncpus=${threads}:mem=800gb

    gunzip -c /media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/data/DEST2_${i}.sync.gz \
    | parallel \
    --jobs ${threads} \
    --pipe \
    -k \
    --cat python3 /media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/scripts/filterSYNCbyName.py \
        --names /media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/data/names_${i}.txt \
        --pools /media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/data/pools_${i}.txt \
        --exclude /media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/data/exclude.txt \
        --input {} \
        --output /media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/data/New_${i} \
       | gzip > /media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/data/DEST2_${i}_sub.sync.gz

    """ >/media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/shell/${i}_SyncFilter.qsub

    gunzip -c /media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/data/DEST2_${i}_sub.sync.gz |
        sort -T/media/inter/mkapun/ -k1,1 -k2,2n \
            >/media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/data/DEST2_${i}_sub.sync

    pigz -f /media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/data/DEST2_${i}_sub.sync

done

## make Names and Pools inputfiles for Grenedalf
awk -F ',' 'NR>1{print $1}' /media/inter/mkapun/projects/ImPoolation/data/meta.csv \
    >/media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/data/names.txt

awk -F ',' 'NR>1{print $24}' /media/inter/mkapun/projects/ImPoolation/data/meta.csv \
    >/media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/data/pools.txt

for i in poolsnp snape; do

    for index in ${!Chrom[@]}; do

        echo """
        #!/bin/sh

        ## name of Job
        #PBS -N ${i}_${Chrom[index]}_grenedalf

        ## Redirect output stream to this file.
        #PBS -o /media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/results/${i}_${Chrom[index]}_log.txt

        ## Stream Standard Output AND Standard Error to outputfile (see above)
        #PBS -j oe

        ## Select a maximum of 100 cores and 800gb of RAM
        #PBS -l select=1:ncpus=100:mem=500gb

        ##### run analyses #####

        # activate the grenedalf environment
        module load Tools/grenedalf-0.1.0

        grenedalf diversity \
            --window-type regions \
            --window-region ${Chrom[index]}":1-"${Lengths[index]} \
            --pool-sizes /media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/data/New_${i}_pools_sub.txt \
            --sample-name-list /media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/data/New_${i}_names_sub.txt \
            --measure all \
            --sync-path /media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/data/DEST2_${i}_sub.sync.gz \
            --out-dir /media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/results/ \
            --file-prefix ${i}"_"${Chrom[index]} \
            --window-region-skip-empty \
            --allow-file-overwriting \
            --threads ${threads} \
            --min-coverage 2 \
            --max-coverage 1000

        """ >/media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/shell/${i}"_"${Chrom[index]}_grenedalf.qsub

        qsub /media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/shell/${i}"_"${Chrom[index]}_grenedalf.qsub

    done

done

## summarize Grenedalf output in tabular form, append info on continent from Metadata file and calculate weighted genome-wide averges
python /media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/scripts/summariseGrenedalf.py \
    --input /media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/results/poolsnp \
    --meta /media/inter/mkapun/projects/DESTv2/Grenedalf_PopGen/data/meta_cov.csv \
    >/media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/results/Grenedalf_poolsnp.summary

## summarize Grenedalf output in tabular form, append info on continent from Metadata file and calculate weighted genome-wide averges
python /media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/scripts/summariseGrenedalf.py \
    --input /media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/results/snape \
    --meta /media/inter/mkapun/projects/DESTv2/Grenedalf_PopGen/data/meta_cov.csv \
    >/media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/results/Grenedalf_snape.summary

## plot Distributions of PopGen stats for continents
echo '''

library(tidyverse)
library(cowplot)

DATA=read.table("/media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/results/Grenedalf_poolsnp.summary",header=T,sep="\t")

DATA.gw.theta_pi.p<-na.omit(DATA) %>%
    filter(DATA$Chrom=="GenomeWide"& DATA$Stat!="snp_count")

P.pool <- ggplot(DATA.gw.theta_pi.p,aes(x=Continent,y=as.numeric(Value)))+
    geom_boxplot()+
    ggtitle("PoolSNP")+
    facet_wrap(.~Stat,scales="free_y")+
    ylab("PopGen Statistic")+
    theme_bw()

DATA=read.table("/media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/results/Grenedalf_snape.summary",header=T,sep="\t")

DATA.gw.theta_pi.s<-na.omit(DATA) %>%
    filter(DATA$Chrom=="GenomeWide"& DATA$Stat!="snp_count")

P.snape <- ggplot(DATA.gw.theta_pi.s,aes(x=Continent,y=as.numeric(Value)))+
    geom_boxplot()+
    ggtitle("SNAPE")+
    facet_wrap(.~Stat,scales="free_y")+
    ylab("PopGen Statistic")+
    theme_bw()

P<-plot_grid(P.pool,P.snape,
    nrow=2)

DATA.gw.theta_pi.s$Caller <- rep("SNAPE",nrow(DATA.gw.theta_pi.s))
DATA.gw.theta_pi.p$Caller <- rep("PoolSNP",nrow(DATA.gw.theta_pi.p))
DATA.new<-spread(rbind(DATA.gw.theta_pi.s,DATA.gw.theta_pi.p),Caller,Value)


P <- ggplot(DATA.new,aes(x=SNAPE,y=PoolSNP,col=Stat))+
    geom_point()+
     facet_wrap(.~Stat,scales="free")+
    geom_smooth(
    method = "lm",
    formula = y ~ x)+
    theme_bw()


ggsave("/media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/results/Grenedalf_corr.pdf",
    P,
    width=16,
    height=8)

ggsave("/media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/results/Grenedalf_corr.png",
    P,
    width=16,
    height=8)


''' >/media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/results/Grenedalf.r

Rscript /media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/results/Grenedalf.r

### Now visualize
echo '''

library(tidyverse)
library(rworldmap)
library(kriging)

for (p in c("poolsnp","snape")){

    DATA=read.table(paste0("/media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/results/Grenedalf_",p,".summary"),header=T,sep="\t")

    DATA$Latitude <- as.factor(round(DATA$Latitude,3))
    DATA$Longitude <- as.factor(round(DATA$Longitude,3))

    for (t in c("theta_pi_abs","theta_watterson_abs","tajimas_d")){
        for (i in c("North_America","Europe","Oceania")){

            DATA.pi.summary<-na.omit(DATA) %>%
                filter(DATA$Chrom=="GenomeWide"& DATA$Stat==t & DATA$Continent==i) %>%
                group_by(Latitude,Longitude) %>%
                summarize(pi=mean(Value),N=n())

            DATA.pi.summary$Latitude <- as.numeric(as.character(DATA.pi.summary$Latitude))
            DATA.pi.summary$Longitude <- as.numeric(as.character(DATA.pi.summary$Longitude))

            X=c(min(DATA.pi.summary$Longitude)-2.5,max(DATA.pi.summary$Longitude)+2.5)
            Y=c(min(DATA.pi.summary$Latitude)-2.5,max(DATA.pi.summary$Latitude)+2.5)
            newmap<-getMap(resolution="low")
            color=colorRampPalette(c("darkgreen","green","white","orange","brown"))

            png(paste0("/media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/results/Grenedalf_",p,"_",i,"_",t,".png"),
                width=1200,
                height=800)

            K=kriging(DATA.pi.summary$Longitude,DATA.pi.summary$Latitude,response=DATA.pi.summary$pi,pixels=250)
            Zi=min(DATA.pi.summary$pi)
            Za=max(DATA.pi.summary$pi)
            par(cex=1.5,mar=c(4,4,2,4))
            image(K,zlim=c(Zi,Za),col=color(100),xlim=X,ylim=Y)
            plot(newmap,col=rgb(0,0,0,0.2),add=T,xlim=X1,ylim=Y1)
            points(DATA.pi.summary$Longitude,DATA.pi.summary$Latitude,pch=16,cex=2)
            #legend.col(color(100),seq(Zi,Za,0.0001))
            dev.off()
        }
    }

}

''' >/media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/results/Grenedalf_visual.r

Rscript /media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen/results/Grenedalf_visual.r
