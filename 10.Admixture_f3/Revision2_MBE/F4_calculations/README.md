Independent estimates of the proportion of African ancestry in northeastern samples using the F4 ratio as implemented in poolfstat of all samples from northeastern America using the following configurations: 

    i) Numerator F4: (FR_Ill_Sai_1_2017-09-16, SIM_SIM_w501_1_NA-MM-DD) ;(ZA_Lim_Pha_1_2010-07-16, Candidate)

    ii) Denominator F4: (FR_Ill_Sai_1_2017-09-16, SIM_SIM_w501_1_NA-MM-DD);(ZA_Lim_Pha_1_2010-07-16,US_Cal_Sto_1_2013-09-03) 

where SIM=outgroup; FR_Ill_Sai_1_2017-09-16=European representative; and US_Cal_Sto_1_2013-09-03="pure" North American representant (note that if this were introduced with African ancestry, we might expect an underestimation of the candidate African ancestry). Also, I chose the African representant as the most pure based on the unsupervised clustering algorithm (à la structure) that MG is developing for PoolSeq data.

The ratio alpha of these 2 F4 allows to estimate the contribution of African ancestry as 1-alpha (see Paterson et. al, 2012 for details or the poolfstat vignette).

