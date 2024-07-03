######
# Setup
#####
library(tidyverse)
library(factoextra)
setwd("C:/Users/David/Desktop/Bergland/data")
df <- read.csv("dest_v2.samps_3May2024.csv")

#####
# Initial data filtering
#####

# Filter `df` to only the samples that passed the quality check to be clustered.
# Part of this quality check is having <5% D. simulans contamination.
df_clust <- df %>% filter(!is.na(df$cluster2.0_k4), 
                          !is.na(df$cluster2.0_k8)) %>%
  group_by(locality) %>%
  slice_max(nFlies, with_ties=FALSE) %>% 
  ungroup() %>%
  select(sampleId, locality, country, nFlies, continent, 
         cluster2.0_k4, cluster2.0_k8) %>%
  mutate(
    cluster_k4_region = case_when(
      cluster2.0_k4 == 1 ~ "Africa",
      cluster2.0_k4 == 2 ~ "Europe_east",
      cluster2.0_k4 == 3 ~ "Europe_west",
      cluster2.0_k4 == 4 ~ "Americas"),
    cluster_k8_region = case_when(
      cluster2.0_k8 == 1 ~ "Africa1",
      cluster2.0_k8 == 2 ~ "Europe_west",
      cluster2.0_k8 == 3 ~ "Africa2",
      cluster2.0_k8 == 4 ~ "Americas",
      cluster2.0_k8 == 5 ~ "Caribbean",
      cluster2.0_k8 == 6 ~ "United_States_east",
      cluster2.0_k8 == 7 ~ "Europe_suture_zone",
      cluster2.0_k8 == 8 ~ "Europe_east"))

#####
# Output of popinfo files
#####

# Output for EU models
df_eu <- df_clust %>% filter(continent == "Europe") %>%
  select(-continent)

# Write to popinfo files
write_tsv(df_eu %>% 
            select(sampleId, cluster_k4_region) %>%
            arrange(cluster_k4_region),
          "k4_Europe.popinfo",
          col_names=FALSE)
write_tsv(df_eu %>% 
            select(sampleId, cluster_k8_region) %>%
            arrange(cluster_k8_region),
          "k8_Europe.popinfo",
          col_names=FALSE)

# Output for Americas models
df_am <- df_clust %>% 
  filter(cluster2.0_k8 %in% c(4, 5, 6)) %>%
  filter(continent != "Oceania") %>%
  select(-continent)

# Write to popinfo file
write_tsv(df_am %>% 
            select(sampleId, cluster_k8_region) %>%
            arrange(cluster_k8_region),
          "k8_Americas.popinfo",
          col_names=FALSE)

# Output for Transatlantic models with clustering at k=4, but Guinean and Zambian
# samples given different, novel clustering IDs so that they can be distinguished
# when parsing the popinfo with moments. Use this popinfo for Transatlantic analyses 
# of whether the American flies result from admixture of (EUE xor EUW) and 
# (Guinea xor Zambia), i.e. four-way model comparison.
df_ta_k4 <- df_clust %>%
  filter(country %in% c("Guinea", "Zambia") |
         continent %in% c("Europe", "North_America", "South_America")) %>%
  mutate(cluster_k4_region = case_when(country == "Guinea" ~ "Guinea",
                                       country == "Zambia" ~ "Zambia",
                                       TRUE ~ cluster_k4_region)) %>%
  select(sampleId, cluster_k4_region) %>%
  arrange(cluster_k4_region)

# Write to popinfo file
write_tsv(df_ta_k4,
          "k4_Transatlantic.popinfo",
          col_names=FALSE)

# Output for Transatlantic models that are the same as the one directly above, but
# "expanded" with additional African samples: Cameroon, Gabon, and Nigeria added
# to Guinea in the Sub-Saharan West, and Zimbabwe added to Zambia
W_Afr_countries <- c("Guinea", "Cameroon", "Gabon", "Nigeria")
E_Afr_countries <- c("Zambia", "Zimbabwe")
df_tae_k4 <- df_clust %>%
  filter(country %in% c(W_Afr_countries, E_Afr_countries) |
         continent %in% c("Europe", "North_America", "South_America")) %>%
  mutate(cluster_k4_region = case_when(country %in% W_Afr_countries ~ "Africa_west",
                                       country %in% E_Afr_countries ~ "Africa_east",
                                       TRUE ~ cluster_k4_region)) %>%
  select(sampleId, cluster_k4_region) %>%
  arrange(cluster_k4_region)

# Write to popinfo file
write_tsv(df_tae_k4,
          "k4_Transatlantic_expandedAfr.popinfo",
          col_names=FALSE)

# Output for "expanded Transatlantic models, but with Australia instead of the
# Americas.
W_Afr_countries <- c("Guinea", "Cameroon", "Gabon", "Nigeria")
E_Afr_countries <- c("Zambia", "Zimbabwe")
df_aus_k4 <- df_clust %>%
  filter(country %in% c(W_Afr_countries, E_Afr_countries) |
           continent %in% c("Europe", "Oceania")) %>%
  mutate(cluster_k4_region = case_when(country %in% W_Afr_countries ~ "Africa_west",
                                       country %in% E_Afr_countries ~ "Africa_east",
                                       continent == "Oceania" ~ "Australia",
                                       TRUE ~ cluster_k4_region)) %>%
  select(sampleId, cluster_k4_region) %>%
  arrange(cluster_k4_region)

# Write to popinfo file
write_tsv(df_aus_k4,
          "k4_Australia.popinfo",
          col_names=FALSE)

#####
# EDA and throwaway code
#####

# GN = Guinea
# ZM = Zambia
# GA = Gabon

# EDA of Africa
df %>% #filter(!is.na(df$cluster2.0_k4)) %>%
  filter(continent == "Africa") %>%
  ggplot(aes(long, lat, color=country)) +
  geom_point()

df %>% filter(country == "Guinea")

# EDA of SE NA
df_clust %>% filter(cluster_k8_region == "United_States_east") %>% 
  pull(locality)
df %>% filter(cluster2.0_k8 == 6) %>% pull(locality)

df %>% filter(cluster2.0_k8 %in% c(4, 5, 6)) %>%
  filter(continent != "Oceania") %>%
  select(locality, lat, long)
ggplot(aes(x=long, y=lat)) +
  geom_point()

df_clust %>% filter(cluster_k8_region == "Europe_suture_zone") %>%
  print(n=45) %>%
  filter(country %in% c("Egypt", "Cyprus", "Turkey")) %>% 
  group_by(country) %>%
  count()

df %>% filter(!is.na(df$cluster2.0_k4), 
              !is.na(df$cluster2.0_k8),
              continent == "Africa") %>%
  ggplot(aes(long, lat, color=(country %in% c("Zambia", "Zimbabwe")))) +
  geom_point()

df %>% filter(!is.na(df$cluster2.0_k4), 
              !is.na(df$cluster2.0_k8),
              continent == "Africa",
              long < 20) %>%
  pull(country)

# Australia
ggplot(data=df %>% filter(continent == "Oceania",
                          !is.na(df$cluster2.0_k4),
                          !is.na(df$cluster2.0_k8)), aes(long, lat)) +
  geom_point()

df %>% filter(continent == "Oceania",
              !is.na(df$cluster2.0_k4),
              !is.na(df$cluster2.0_k8)) %>%
  select(locality, nFlies)

ggplot(data=df %>% filter(!is.na(df$cluster2.0_k4),
                          !is.na(df$cluster2.0_k8)), aes(long, lat, color=continent)) +
         geom_point()
       