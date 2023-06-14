
library(data.table)


#### Contaminantion --->
contam<-get(load("/yey2sn/DEST2_analysis/filtering/Mean.contamination.Final.Rdata"))
#### DuplicateRates --->
duprate<-get(load("/yey2sn/DEST2_analysis/filtering/DuplicateRates.all.Rdata"))
#### PnPs ---->
#system("wget https://raw.githubusercontent.com/DEST-bio/DESTv2_data_paper/main/misc/PnPs/results/classify_pops.txt?token=GHSAT0AAAAAAB7U52CNPDJUDB7GWR4N27GQZCAF6WA -O pnps.predict.txt")
pnps <-fread("pnps.predict.txt")
######
