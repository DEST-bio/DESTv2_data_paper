library(tidyverse)
library(knitr)

meta <- read.csv("/media/inter/mkapun/projects/DESTv2_data_paper/16.Inversions/data/meta.csv",
    header = T
)

# POOL-SNP
nhm_inversion <- read.table("/media/inter/mkapun/projects/DESTv2_data_paper/16.Inversions/results/frequencies/DEST2_inversions.af",
    header = T,
    na.string = "NA"
)

data <- merge(meta, nhm_inversion, by = "sampleId")
data$continent[data$continent == "Asia"] <- "Europe"
data$continent <- as.factor(data$continent)

# Inversions
inversions <- c("In.3R.Payne", "In.2L.t", "In.2R.Ns", "In.3R.C", "In.3R.K", "In.3R.Mo", "In.3L.P")


dir.create("/media/inter/mkapun/projects/DESTv2_data_paper/16.Inversions/results/clines/")

sink("/media/inter/mkapun/projects/DESTv2_data_paper/16.Inversions/results/clines/InversionStats.txt")
for (inv in inversions) {
    for (cont in levels(data$continent)) {
        newdat <- data %>%
            filter(continent == cont)

        A <- anova(lm(asin(sqrt(newdat[[inv]])) ~ newdat$lat * newdat$lon * as.factor(newdat$year))) %>%
            kable(caption = paste0(inv, "_", cont))
        print(A)
    }
}
sink()
