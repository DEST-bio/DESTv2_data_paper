# Authors: M Kapun & S Steindl
# Analysis of Inversion Markers of the DESTv2 data

# data availability: http://berglandlab.uvadcos.io/gds/

### get data
WD=/media/inter/mkapun/projects/DESTv2_data_paper

mkdir ${WD}/16.Inversions/data
cd ${WD}/16.Inversions/data
wget -O DEST.vcf.gz http://berglandlab.uvadcos.io/vcf/dest.all.PoolSNP.001.50.3May2024.ann.vcf.gz
wget -O meta.csv https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_3May2024.csv

### Convert VCF to SYNC
gunzip -c /media/inter/mkapun/projects/InvChapter/data/DEST.vcf.gz |
    parallel \
        --jobs 80 \
        --pipe \
        -k \
        --cat python3 ${WD}/16.Inversions/scripts/VCF2sync.py \
        --input {} |
    gzip >/media/inter/mkapun/projects/InvChapter/data/DEST.sync.gz

mkdir -p ${WD}/16.Inversions/results/frequencies

### Get count at inversion specific marker SNP
gunzip -c /media/inter/mkapun/projects/InvChapter/data/DEST.sync.gz |
    parallel \
        --pipe \
        --jobs 200 \
        -k \
        --cat python3 ${WD}/16.Inversions/scripts/overlap_in_SNPs.py \
        --source ${WD}/16.Inversions/data/inversion_markers_v6.txt_pos \
        --target {} \
        >${WD}/16.Inversions/results/frequencies/inversions.sync

### Get names of samples in correct order
NAMES=$(gunzip -c /media/inter/mkapun/projects/InvChapter/data/DEST.vcf.gz | head -500 | awk '/^#C/' | cut -f10- | tr '\t' ',')

# Calculate average frequencies for marker SNPs
python3 ${WD}/16.Inversions/scripts/inversion-freqs.py \
    ${WD}/16.Inversions/data/inversion_markers_v6.txt \
    ${WD}/16.Inversions/results/frequencies/inversions.sync $NAMES \
    >${WD}/16.Inversions/results/frequencies/DEST2_inversions.af

echo """
library(tidyverse)
library(knitr)

meta <- read.csv('${WD}/16.Inversions/data/meta.csv',
    header = T
)

# POOL-SNP
nhm_inversion <- read.table('${WD}/16.Inversions/results/frequencies/DEST2_inversions.af',
    header = T,
    na.string = 'NA'
)

data <- merge(meta, nhm_inversion, by = 'sampleId')
data\$continent[data\$continent == 'Asia'] <- 'Europe'
data\$continent <- as.factor(data\$continent)

# Inversions
inversions <- c('In.3R.Payne', 'In.2L.t', 'In.2R.Ns', 'In.3R.C', 'In.3R.K', 'In.3R.Mo', 'In.3L.P')

dir.create('${WD}/16.Inversions/results/clines/')

sink('${WD}/16.Inversions/results/clines/InversionStats.txt')
for (inv in inversions) {
    for (cont in levels(data\$continent)) {
        newdat <- data %>%
            filter(continent == cont)

        A <-anova(lm(asin(sqrt(newdat[[inv]])) ~ newdat\$lat * newdat\$lon * as.factor(newdat\$year))) %>%
            kable(caption = paste0(inv, '_', cont))
        print(A)
    }
}
sink()

""" >${WD}/16.Inversions/scripts/InversionStats.R

Rscript ${WD}/16.Inversions/scripts/InversionStats.R


### now combine with previous estimates
cd ${WD}/16.Inversions/data
## download 'https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fmec.14871&file=mec14871-sup-0001-TableS1.xlsx' as KapunFlatt.xlsx

echo """

library(tidyverse)
library(knitr)
library(readxl)

setwd(${WD})

meta <- read.csv('16.Inversions/data/meta.csv',
    header = T
)

# POOL-SNP
nhm_inversion <- read.table('16.Inversions/results/frequencies/DEST2_inversions.af',
    header = T,
    na.string = 'NA'
)

data <- merge(meta, nhm_inversion, by = 'sampleId') %>%
    select(lat,long,`In.2L.t`,`In.2R.Ns`,`In.3L.P`,`In.3R.Payne`)
colnames(data)<-c('Latitude','Longitude','In(2L)t','In(2R)Ns','In(3L)P','In(3R)Payne')
data\$DEST<-rep(22,nrow(data))

KapFlatt <- read_excel('16.Inversions/data/KapunFlatt.xlsx', skip = 5) %>%
    filter(!(grepl('Kapun',Reference) | grepl('DrosEU',Reference))) %>%
    select(`Lat(N)`,`Long(E)`,`IN(2L)t`,`IN(2R)NS`,`IN(3L)P`,`IN(3R)P`)
colnames(KapFlatt)<-c('Latitude','Longitude','In(2L)t','In(2R)Ns','In(3L)P','In(3R)Payne')
KapFlatt\$DEST<-rep(21,nrow(KapFlatt))
NewDat<-rbind(data,KapFlatt)

write.table(file='FIGUREs/FIGURE2_Invs/FullInvDestv2.txt',
    NewDat,
    quote=FALSE,
    row.names=FALSE)
""" > ${WD}/16.Inversions/scripts/InverionMerge.R

Rscript ${WD}/16.Inversions/scripts/InversionMerge.R

### now plot

Rscript ${WD}/FIGUREs/FIGURE2_Invs/reproduce_fig2_PlotDots.r