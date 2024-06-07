library(tidyverse)
library(knitr)
library(readxl)

setwd("/media/inter/mkapun/projects/DESTv2_data_paper")

meta <- read.csv("16.Inversions/data/meta.csv",
    header = T
)

# POOL-SNP
nhm_inversion <- read.table("16.Inversions/results/frequencies/DEST2_inversions.af",
    header = T,
    na.string = "NA"
)

data <- merge(meta, nhm_inversion, by = "sampleId") %>%
    select(lat, long, `In.2L.t`, `In.2R.Ns`, `In.3L.P`, `In.3R.Payne`)
colnames(data) <- c("Latitude", "Longitude", "In(2L)t", "In(2R)Ns", "In(3L)P", "In(3R)Payne")
data$DEST <- rep(22, nrow(data))

KapFlatt <- read_excel("16.Inversions/data/KapunFlatt.xlsx", skip = 5) %>%
    filter(!(grepl("Kapun", Reference) | grepl("DrosEU", Reference))) %>%
    select(`Lat(N)`, `Long(E)`, `IN(2L)t`, `IN(2R)NS`, `IN(3L)P`, `IN(3R)P`)
colnames(KapFlatt) <- c("Latitude", "Longitude", "In(2L)t", "In(2R)Ns", "In(3L)P", "In(3R)Payne")
KapFlatt$DEST <- rep(21, nrow(KapFlatt))
NewDat <- rbind(data, KapFlatt)

write.table(
    file = "FIGUREs/FIGURE2_Invs/FullInvDestv2.txt",
    NewDat,
    quote = FALSE,
    row.names = FALSE
)
