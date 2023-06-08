#### ---> Make collapse recomendations
#### 

samps <- fread("../FIGURE1/code/dest_v2.samps_25Feb2023.csv")
QC.recs <- fread("QC.recomendations.csv")

full_join(samps,QC.recs) %>%
  mutate(merger_id = paste(city, loc_rep, sep = "_" )) ->
  samps.QC.merger