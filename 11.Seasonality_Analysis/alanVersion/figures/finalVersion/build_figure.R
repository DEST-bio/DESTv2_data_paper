### library
  library(data.table)
  library(ggplot2)
  library(patchwork)

### load data
  ### A, B, C: this file loads the output of the WZA on the un-permuted data. WZA for XtXst pvalue (Chisq based), Contrast pvalue, and GLM pvalue
    load(file="~/XtX_C2_glm.windows.Rdata") ### `win.out`


  ### B: Window permutations?

  ### C: this file contains the GLM permutations used to draw significance threshold for panel C.
    load(file="~/destv2_seasonality_perm.Rdata")) ### `win.all.ag2`

  ###
