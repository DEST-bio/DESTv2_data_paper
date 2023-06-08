#### ---> Make collapse recomendations
#### 

samps <- fread("../FIGURE1/code/dest_v2.samps_25Feb2023.csv")
QC.recs <- fread("QC.recomendations.csv")

full_join(samps,QC.recs) %>%
  mutate(merger_id = paste(city, loc_rep, sep = "_" )) %>%
  filter(collector == "Fournier-Level et al")->
  samps.QC.merger.Fournier

save(samps.QC.merger.Fournier, file = "samps.QC.merger.Fournier.Rdata")

samps.QC.merger.Fournier %>%
  dplyr::select(sampleId,merger_id,Recomendation  ) ->
  master_merger_Fournier


master_merger_Fournier %>%
  group_by(merger_id) %>%
  filter(Recomendation != "High Contamination" ) %>%
  summarise(N = n())


