#### Reconstruct Figure 7 as. of Feb 26, 2024
#### Jcbn

library(tidyverse)
library(reshape2)

#XtX_dat = get(load("/netfiles/nunezlab/Drosophila_resources/Datasets/2023.DEST.2.0._release/Seasonality_Analyses/dest2_glm_baypass_annotation.Rdata"))

summary_win = get(load("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/DESTv2_data_paper/FIGUREs/FIGURE7_seasonal/data_for_reproduction/XtX_C2_glm.windows.Rdata"))


### Make wZA panels
summary_win %>%
  dplyr::select(chr, pos_mean, xtx.wZa.p, C2.wZa.p, wZa.p) %>%
  melt(id=c("pos_mean", "chr")) %>%
  ggplot(aes(x=pos_mean/1e6,y=-value)) +
  geom_line() +
  facet_grid(variable~chr, scales = "free") +
  theme_bw() ->
  all.wZa.p.plot

ggsave(all.wZa.p.plot, file = "all.wZa.p.plot.pdf",
       h = 6.0, w = 10)


