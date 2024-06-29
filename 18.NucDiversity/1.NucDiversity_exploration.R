####
library(tidyverse)
library(data.table)

load("Figure1.RData")
samps <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_8Jun2023.csv")

names(pop_genome_stats_pi_m)[1] = "sampleId"

pop_genome_stats_pi_m %>% 
  filter(variable == "npStat") %>% 
  left_join(samps) ->
  pop_genome_stats_pi_m.samps

pop_genome_stats_pi_m.samps %>%
  filter(lat > 23.5 & lat < 66.5) %>%
  ggplot(
    aes(
      x=lat,
      y=value,
      color=continent
    )
  ) + geom_point() +
  ylab(expression(pi)) +
  geom_smooth(method = "lm", color = "black") ->
  pi_lat
ggsave(pi_lat, file = "pi_lat.pdf")

pop_genome_stats_pi_m.samps %>%
  filter(lat > 23.5 & lat < 66.5) %>%
  ggplot(
    aes(
      x=long,
      y=value,
      color=continent
    )
  ) + geom_point() +
  ylab(expression(pi)) +
  geom_smooth(method = "lm", color = "black") ->
  pi_long
ggsave(pi_long, file = "pi_long.pdf")

pop_genome_stats_pi_m.samps %>%
  filter(lat > 23.5 & lat < 66.5) ->
  temp_pi_samps

cor.test(~ value + lat, data = temp_pi_samps)
cor.test(~ value + long, data = temp_pi_samps)

cor.test(~ value + lat, data = filter(temp_pi_samps, cluster2.0_k4 == 2))
cor.test(~ value + lat, data = filter(temp_pi_samps, cluster2.0_k4 == 3))

cor.test(~ value + long, data = filter(temp_pi_samps, cluster2.0_k4 == 2))
cor.test(~ value + long, data = filter(temp_pi_samps, cluster2.0_k4 == 3))

cor.test(~ value + lat, data = filter(temp_pi_samps, cluster2.0_k4 == 4))
cor.test(~ value + long, data = filter(temp_pi_samps, cluster2.0_k4 == 4))


pop_genome_stats_pi_m.samps %>%
  filter(lat < -23.5 & lat > -66.5) ->
  temp_pi_samps_S

cor.test(~ value + lat, data = filter(temp_pi_samps_S, continent == "Oceania"))
cor.test(~ value + long, data = filter(temp_pi_samps_S, continent == "Oceania"))

cor.test(~ value + lat, data = filter(temp_pi_samps_S, continent == "South_America"))
cor.test(~ value + long, data = filter(temp_pi_samps_S, continent == "Oceania"))


