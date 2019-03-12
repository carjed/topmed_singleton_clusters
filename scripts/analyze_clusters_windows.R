#-----------------------------------------------------------------------------
# NOT RUN: Run per genomic window instead of per sample
#-----------------------------------------------------------------------------
bin_counts <- idsites %>% 
  group_by(BIN) %>% 
  summarise(tot=n())

run_sites <- merge(idsites, bin_counts, by="BIN") %>%
  dplyr::filter(tot>1000) %>%
  dplyr::filter(D2>0 & D2<2e6) %>%
  ungroup()

bins_part <- multidplyr::partition(run_sites, BIN, cluster=cluster)
bins_part %>%
  cluster_library("tidyverse") %>%
  cluster_library("mixtools") %>%
  cluster_library("MASS") %>%
  cluster_assign_value("getExpLogLik", getExpLogLik) %>%
  cluster_assign_value("getNR2", getNR2) %>%
  cluster_assign_value("fitExpMix", fitExpMix) %>%
  cluster_assign_value("runExpMix", runExpMix)

exp_fits_b <- bins_part %>%
  do(fitExpMix(x=.$D2, scale=10, mincomp=3, maxcomp=4)) %>%
  as_tibble()