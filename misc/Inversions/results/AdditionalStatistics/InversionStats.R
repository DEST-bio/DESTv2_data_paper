library(tidyverse)
library(knitr)

meta <- read.csv("/media/inter/mkapun/projects/DESTv2_data_paper/misc/Inversions/data/meta.csv",
    header = T
)

# POOL-SNP
nhm_inversion <- read.table("/media/inter/mkapun/projects/DESTv2_data_paper/misc/Inversions/results/PoolSNP_nhm_inversion.af", header = T, na.string = "na")

data <- merge(meta, nhm_inversion, by = "Sample")
data$continent[data$continent == "Asia"] <- "Europe"
data$continent <- as.factor(data$continent)

# Inversions
inversions <- c("In.3R.Payne", "In.2L.t", "In.2R.Ns", "In.3R.C", "In.3R.K", "In.3R.Mo", "In.3L.P")


sink("/media/inter/mkapun/projects/DESTv2_data_paper/misc/Inversions/results/AdditionalStatistics/InversionStats.txt")
for (inv in inversions) {
    for (cont in levels(data$continent)) {
        newdat <- data %>%
            filter(continent == cont)

        A=anova(lm(newdat[[inv]] ~ newdat$lat * newdat$lon)) %>%
            kable(caption = paste0(inv, "_", cont)))
    }
}
sink()
