## get SYNC data

PWD=/media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen

mkdir ${PWD}/data

## Convert PoolSNP VCF 2 SYNC
gunzip -c /media/inter/ssteindl/DEST/DEST2_NHM/collapsed/PoolSNP/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.vcf.gz |
    parallel \
        --jobs 200 \
        --pipe \
        -k \
        --cat python3 ${PWD}/scripts/VCF2sync.py \
        --input {} |
    gzip >${PWD}/data/DEST2_poolsnp.sync.gz

## Copy meta data
cp /media/inter/mkapun/projects/ImPoolation/data/meta.csv ${PWD}/data

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
ZW_Mat_Vic_1_2001-07-16''' >${PWD}/data/exclude.txt

## Define Chromsomal arms to consider
Chrom=("4" "2L" "2R" "3L" "3R" "X")
Lengths=(1348131 23513712 25286936 28110227 32079331 23542271)
threads=50

## make Names and Pools inputfiles for Grenedalf
awk -F ',' 'NR>1{print $1}' ${PWD}/data/meta.csv \
    >${PWD}/data/names.txt

awk -F ',' 'NR>1{print $24}' ${PWD}/data/meta.csv \
    >${PWD}/data/pools.txt

## make subset of Sync based on populations to exclude
echo """
    #!/bin/sh

    ## name of Job
    #PBS -N mapping

    ## Redirect output stream to this file.
    #PBS -o ${PWD}/results/poolsnp_filterSync_log.txt

    ## Stream Standard Output AND Standard Error to outputfile (see above)
    #PBS -j oe

    ## Select a maximum of ${threads} cores and 800gb of RAM
    #PBS -l select=1:ncpus=${threads}:mem=800gb

    gunzip -c ${PWD}/data/DEST2_poolsnp.sync.gz \
    | parallel \
    --jobs ${threads} \
    --pipe \
    -k \
    --cat python3 ${PWD}/scripts/filterSYNCbyName.py \
        --names ${PWD}/data/names.txt \
        --pools ${PWD}/data/pools.txt \
        --exclude ${PWD}/data/exclude.txt \
        --input {} \
        --output ${PWD}/data/New_poolsnp \
       | gzip > ${PWD}/data/DEST2_poolsnp_sub.sync.gz

    """ >${PWD}/shell/poolsnp_SyncFilter.qsub

sh ${PWD}/shell/poolsnp_SyncFilter.qsub

##
gunzip -c ${PWD}/data/DEST2_poolsnp_sub.sync.gz |
    sort -T/media/inter/mkapun/ -k1,1 -k2,2n |
    pigz -f >${PWD}/data/DEST2_poolsnp_sub.sync.gz

## make names input for Grenedalf
rm -rf ${PWD}/data/poolsnp_names.txt
var=1
while IFS=$"" read -r column1; do
    printf "DEST2_poolsnp_sub."$var"\t"$column1"\n" >>${PWD}/data/poolsnp_names.txt
    ((var++))
done <${PWD}/data/New_poolsnp_names_sub.txt

## make pools input for Grenedalf
rm -rf ${PWD}/data/poolsnp_pools.txt
paste ${PWD}/data/New_poolsnp_names_sub.txt \
    ${PWD}/data/New_poolsnp_pools_sub.txt \
    >>${PWD}/data/poolsnp_pools.txt

## run GRENEDALF per chromosomes
echo """
#!/bin/sh

## name of Job
#PBS -N poolsnp_${Chrom[index]}_grenedalf

## Redirect output stream to this file.
#PBS -o ${PWD}/results/poolsnp_${Chrom[index]}_log.txt

## Stream Standard Output AND Standard Error to outputfile (see above)
#PBS -j oe

## Select a maximum of 100 cores and 800gb of RAM
#PBS -l select=1:ncpus=100:mem=500gb

##### run analyses #####

# activate the grenedalf environment
module load Tools/grenedalf-0.3.0

grenedalf diversity \
    --window-type chromosomes \
    --pool-sizes ${PWD}/data/poolsnp_pools.txt \
    --rename-samples-file ${PWD}/data/poolsnp_names.txt \
    --measure all \
    --sync-path ${PWD}/data/DEST2_poolsnp_sub.sync.gz \
    --out-dir ${PWD}/results/ \
    --file-prefix poolsnp_chrom \
    --allow-file-overwriting \
    --threads ${threads} \
    --filter-sample-min-coverage 2 \
    --filter-sample-max-coverage 250 \
    --filter-sample-min-count 2 

""" >${PWD}/shell/poolsnp_chrom_grenedalf.qsub

sh ${PWD}/shell/poolsnp_chrom_grenedalf.qsub

## run GRENEDALF per window
for window in 10000 50000 100000; do

    echo """
    #!/bin/sh

    ## name of Job
    #PBS -N poolsnp_${Chrom[index]}_grenedalf

    ## Redirect output stream to this file.
    #PBS -o ${PWD}/results/poolsnp_${Chrom[index]}_log.txt

    ## Stream Standard Output AND Standard Error to outputfile (see above)
    #PBS -j oe

    ## Select a maximum of 100 cores and 800gb of RAM
    #PBS -l select=1:ncpus=100:mem=500gb

    ##### run analyses #####

    # activate the grenedalf environment
    module load Tools/grenedalf-0.3.0

    grenedalf diversity \
        --window-type sliding \
        --window-sliding-width ${window} \
        --pool-sizes ${PWD}/data/poolsnp_pools.txt \
        --rename-samples-file ${PWD}/data/poolsnp_names.txt \
        --measure all \
        --sync-path ${PWD}/data/DEST2_poolsnp_sub.sync.gz \
        --out-dir ${PWD}/results/ \
        --file-prefix poolsnp_${window} \
        --allow-file-overwriting \
        --threads ${threads} \
        --filter-sample-min-coverage 2 \
        --filter-sample-max-coverage 250 \
        --filter-sample-min-count 2 

    """ >${PWD}/shell/poolsnp_${window}_grenedalf.qsub

    sh ${PWD}/shell/poolsnp_${window}_grenedalf.qsub &

done

## compress
pigz -f ${PWD}/results/poolsnp_*diversity.csv

## get URLs for BAM files
cd ${PWD}/data
curl -O http://berglandlab.uvadcos.io/dest_mapped/local.addresses.DEST2.0.txt

## download all BED files with positions to exclude
mkdir ${PWD}/data/BED
cd ${PWD}/data/BED
while IFS=$"" read -r ID; do
    if [[ "$ID" == *".mel.bam"* ]]; then
        curl -O ${ID%%.mel.bam*}.bed.gz
    fi
done <${PWD}/data/local.addresses.DEST2.0.txt

## calculate windows to be excluded for each BED file and window size
mkdir ${PWD}/data/BED/windows

for FILE in ${PWD}/data/BED/*.bed.gz; do
    for WINDOW in 10000 50000 100000; do

        tmp=${FILE##*/}
        ID=${tmp%%.bed.gz*}

        echo """
        #!/bin/sh

        ## name of Job
        #PBS -N Windows

        ## Redirect output stream to this file.
        #PBS -o ${PWD}/data/logs

        ## Stream Standard Output AND Standard Error to outputfile (see above)
        #PBS -j oe

        ## Select a maximum of 100 cores and 800gb of RAM
        #PBS -l select=1:ncpus=1:mem=5gb

        ##### run analyses #####

        # activate the grenedalf environment
        module load Tools/grenedalf-0.3.0
        python ${PWD}/scripts/BED2Window.py \
            --window ${WINDOW} \
            --input ${FILE} \
            > ${PWD}/data/BED/windows/${ID}_${WINDOW}.txt

        """ >${PWD}/data/logs/${ID}_${WINDOW}.qsub
        qsub ${PWD}/data/logs/${ID}_${WINDOW}.qsub
    done

done

## concatenate files across all samples
for WINDOW in 10000 50000 100000; do
    cat ${PWD}/data/BED/windows/*_${WINDOW}.txt |
        gzip >>${PWD}/data/Window_${WINDOW}.txt.gz
done

for FILE in ${PWD}/data/BED/*.bed.gz; do

    tmp=${FILE##*/}
    ID=${tmp%%.bed.gz*}

    echo """
    #!/bin/sh

    ## name of Job
    #PBS -N Windows

    ## Redirect output stream to this file.
    #PBS -o ${PWD}/data/logs

    ## Stream Standard Output AND Standard Error to outputfile (see above)
    #PBS -j oe

    ## Select a maximum of 100 cores and 800gb of RAM
    #PBS -l select=1:ncpus=1:mem=2gb

    ##### run analyses #####

    # activate the grenedalf environment
    module load Tools/grenedalf-0.3.0
    python ${PWD}/scripts/BED2Window.py \
        --window chrom \
        --input ${FILE} \
        > ${PWD}/data/BED/windows/${ID}_chrom.txt

    """ >${PWD}/data/logs/${ID}_chrom.qsub
    qsub ${PWD}/data/logs/${ID}_chrom.qsub

done

## concatenate files across all samples
cat ${PWD}/data/BED/windows/*_chrom.txt |
    gzip >>${PWD}/data/Window_chrom.txt.gz

# Summarize Grendalf for each chromosome
printf "ID\tContinent\tCountry\tLatitude\tLongitude\tChrom\tPos\tStat\tValue\n" |
    gzip >${PWD}/results/Grenedalf_poolsnp_chrom.summary.gz

python ${PWD}/scripts/summariseGrenedalf.py \
    --input ${PWD}/results/poolsnp_chromdiversity.csv.gz \
    --coverage ${PWD}/data/Window_chrom.txt.gz \
    --meta ${PWD}/data/meta.csv \
    --chrom |
    gzip >>${PWD}/results/Grenedalf_poolsnp_chrom.summary.gz

# Summarize Grendalf in windows (10k, 50k and 100k)
for i in 10000 50000 100000; do

    printf "ID\tContinent\tCountry\tLatitude\tLongitude\tChrom\tPos\tStat\tValue\n" |
        gzip >${PWD}/results/Grenedalf_poolsnp_${i}.summary.gz

    python ${PWD}/scripts/summariseGrenedalf.py \
        --input ${PWD}/results/poolsnp_${i}diversity.csv.gz \
        --coverage ${PWD}/data/Window_${i}.txt.gz \
        --meta ${PWD}/data/meta.csv |
        gzip >>${PWD}/results/Grenedalf_poolsnp_${i}.summary.gz
done

## plot Distributions of PopGen stats for continents
echo """

library(tidyverse)
library(cowplot)

setwd('${PWD}')
#setwd('/media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen')

## read Chromosome-wide data
DATA=read.table('results/Grenedalf_poolsnp_chrom.summary.gz',header=T,sep="\t")

## based on the final Metadata file from the MS (Supp. Table S1)
Filter=read.table('data/filterIDs.txt',header=T,sep="\t")
Keep <- Filter %>%
    filter(Recommendation == 'Pass')

DATA.genom <-DATA %>%
    group_by(Continent,Latitude, Longitude,Stat,ID) %>%
    summarize(Value=mean(Value)) %>%
    filter(ID %in% Keep$sampleId)

P.pool <- ggplot(DATA.genom,
    aes(x=Continent,
        y=as.numeric(Value)))+
    geom_boxplot()+
    ggtitle('PoolSNP')+
    facet_wrap(.~Stat,
        scales='free_y')+
    ylab('PopGen Statistic')+
    theme_bw()
P.pool

# DATA.gw.theta_pi.s$Caller <- rep("SNAPE",nrow(DATA.gw.theta_pi.s))
# DATA.gw.theta_pi.p$Caller <- rep("PoolSNP",nrow(DATA.gw.theta_pi.p))
# DATA.new<-spread(rbind(DATA.gw.theta_pi.s,DATA.gw.theta_pi.p),Caller,Value)


# P <- ggplot(DATA.new,aes(x=SNAPE,y=PoolSNP,col=Stat))+
#     geom_point()+
#      facet_wrap(.~Stat,scales="free")+
#     geom_smooth(
#     method = "lm",
#     formula = y ~ x)+
#     theme_bw()

dir.create('results/figures')
ggsave('results/figures/Grenedalf_Stats.pdf',
    P.pool,
    width=16,
    height=8)

ggsave('results/figures/Grenedalf_Stats.png',
    P.pool,
    width=16,
    height=8)

""" >${PWD}/results/Grenedalf.r

Rscript ${PWD}/results/Grenedalf.r

### Now visualize by kriging
# echo """

# library(tidyverse)
# library(rworldmap)
# library(kriging)

# setwd('${PWD}')
# #setwd('/media/inter/mkapun/projects/DESTv2_data_paper/misc/Grenedalf_PopGen')

# ## read Chromosome-wide data
# DATA=read.table('results/Grenedalf_poolsnp_chrom.summary.gz',header=T,sep="\t")

# ## based on the final Metadata file from the MS (Supp. Table S1)
# Filter=read.table('data/filterIDs.txt',header=T,sep="\t")
# Keep <- Filter %>%
#     filter(Recommendation == 'Pass')

# DATA.genom <-DATA %>%
#     group_by(Continent,Latitude, Longitude,Stat,ID) %>%
#     summarize(Value=mean(Value)) %>%
#     filter(ID %in% Keep$sampleId)

# DATA.genom$Latitude <- as.factor(round(DATA.genom$Latitude,3))
# DATA.genom$Longitude <- as.factor(round(DATA.genom$Longitude,3))

# for (t in c('theta_pi_rel','theta_watterson_rel','tajimas_d')){
#     for (i in c('North_America','Europe','Oceania')){

#         DATA.pi.summary<-DATA.genom %>%
#             filter(Stat==t & Continent==i) %>%
#             group_by(Latitude,Longitude) %>%
#             summarize(pi=mean(Value),N=n())

#         DATA.pi.summary$Latitude <- as.numeric(as.character(DATA.pi.summary$Latitude))
#         DATA.pi.summary$Longitude <- as.numeric(as.character(DATA.pi.summary$Longitude))

#         X=c(min(DATA.pi.summary$Longitude)-2.5,max(DATA.pi.summary$Longitude)+2.5)
#         Y=c(min(DATA.pi.summary$Latitude)-2.5,max(DATA.pi.summary$Latitude)+2.5)
#         newmap<-getMap(resolution='low')
#         color=colorRampPalette(c('darkgreen','green','white','orange','brown'))

#         png(paste0('results/figures_',i,'_',t,'.png'),
#             width=1200,
#             height=800)

#         K=kriging(DATA.pi.summary$Longitude,DATA.pi.summary$Latitude,response=DATA.pi.summary$pi,pixels=250)
#         Zi=min(DATA.pi.summary$pi)
#         Za=max(DATA.pi.summary$pi)
#         par(cex=1.5,mar=c(4,4,2,4))
#         image(K,zlim=c(Zi,Za),col=color(100),xlim=X,ylim=Y)
#         plot(newmap,col=rgb(0,0,0,0.2),add=T,xlim=X1,ylim=Y1)
#         points(DATA.pi.summary$Longitude,DATA.pi.summary$Latitude,pch=16,cex=2)
#         #legend.col(color(100),seq(Zi,Za,0.0001))
#         dev.off()
#         }
#     }

# """ >${PWD}/results/Grenedalf_visual.r

# Rscript ${PWD}/results/Grenedalf_visual.r
