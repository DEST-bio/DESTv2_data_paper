fl_shufAll <- list.files("/scratch/aob2x/dest2_baypass/contrast_perms/", "contrast_perm_subpool_", full.name=T) ### this was totally random shuffling
fl_locYear <- list.files("/scratch/aob2x/dest2_baypass/contrast_perms_v3", "contrast_perm_subpool_", full.name=T); fl_locYear<-fl_locYear[!grepl("RAW", fl_locYear)] ### this shuffles within locality_year
