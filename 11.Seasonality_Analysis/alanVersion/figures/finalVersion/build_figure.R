### library
  library(data.table)
  library(ggplot2)
  library(patchwork)

### load data
  setwd("~/")

  ### A, B, C: this file loads the output of the WZA on the un-permuted data. WZA for XtXst pvalue (Chisq based), Contrast pvalue, and GLM pvalue
    load(file="XtX_C2_glm.windows.Rdata") ### `win.out`

  ### B: C2 permutations from WZA to build threshold?

  ### C: this file contains the GLM permutations used to draw significance threshold for panel C.
    load(file="destv2_seasonality_perm.Rdata")) ### `win.all.ag2`

  ### D: this file contains the seasonal pairs. reference this file: `DESTv2_data_paper/11.Seasonality_Analysis/alanVersion/figures/mapFigure.R`
    load("seasonalpair.pca.meta.Rdata")


  ### E: this file should contain what we want
    load("enrichment.NoCore20_seas_europe_yearPop_Ran.Rdata")

  ### F: this file contains the detailed SNP level metrics that allow us to build the hexbin plot. Also contains the thresholds for the POD analysis of C2 and XTX
    load(file="dest2_glm_baypass_annotation_pod.Rdata")
