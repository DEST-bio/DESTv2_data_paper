# Seasonality Visualization

```bash
/netfiles/nunezlab/Drosophila_resources/Datasets/2023.DEST.2.0._release/Seasonality_Analyses/dest2_glm_baypass_annotation.Rdata

/netfiles/nunezlab/Drosophila_resources/Datasets/2023.DEST.2.0._release/Seasonality_Analyses/XtX_C2_glm.windows.Rdata
```
load objects
```r
library(tidyverse)
library(reshape2)

XtX_dat = get(load("/netfiles/nunezlab/Drosophila_resources/Datasets/2023.DEST.2.0._release/Seasonality_Analyses/dest2_glm_baypass_annotation.Rdata"))

summary_win = get(load("/netfiles/nunezlab/Drosophila_resources/Datasets/2023.DEST.2.0._release/Seasonality_Analyses/XtX_C2_glm.windows.Rdata"))
```
create compound object for the wZA stats.
```r
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

```


explore objects
```r
ggplot() +
geom_line(
data=summary_win,
aes(x=pos_mean/1e6,y=-xtx.wZa.p)) +
facet_grid(~chr, scales = "free_x") +
theme_bw() ->
xtx.wZa.p.plot

ggsave(xtx.wZa.p.plot, file = "xtx.wZa.p.plot.pdf",
h = 3, w = 12)

summary_win %>% group_by(chr) %>% slice_min(xtx.wZa.p, n=1) %>% as.data.frame %>% select(chr,pos_mean, xtx.genes, xtx.wZa.p) -> tops.xtx
```

```r
ggplot() +
geom_line(
data=summary_win,
aes(x=pos_mean/1e6,y=-C2.wZa.p)) +
facet_grid(~chr, scales = "free_x") +
theme_bw() ->
C2.wZa.p.plot

ggsave(C2.wZa.p.plot, file = "C2.wZa.p.plot.pdf",
h = 3, w = 12)

summary_win %>% group_by(chr,invName) %>% slice_min(C2.wZa.p, n=3) %>% as.data.frame %>% select(chr,pos_mean,invName, c2.genes, C2.wZa.p) -> tops.c

tops.c %>%
filter(chr == "2L")


tops.c %>%
filter(chr == "3L")

```