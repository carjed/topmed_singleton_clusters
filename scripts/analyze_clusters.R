##############################################################################
# load libraries
##############################################################################
library(tidyverse)
library(scales)
library(MASS)
library(mixtools)
library(broom)
library(ggExtra)
library(gganimate)
library(depmixS4)
library(forecast)
library(arrangements)

# devtools::install_github("hadley/multidplyr")
library(multidplyr)
library(viridis)

options(dplyr.width = Inf)

source("/mnt/norbert/home/jedidiah/projects/ERV_mutation_hotspots/scripts/cluster_functions.R")

##############################################################################
# Read and prep data
##############################################################################

# African, CSAsian, Easian, European, Native American, Oceania, Middle Eastern
anc_header <- c("ID", "AFR", "CSA", "EAS", "EUR", "NAT", "OC", "ME")

# deprecated build37 ancestry for freeze 3 data
# anc <- read_delim("/mnt/norbert/data/topmed/ancestry/global_list_19815.txt", delim=" ", col_names=F) 

# freeze 5 ancestry
anc <- read_delim("/mnt/norbert/data/topmed/ancestry/global_list_b38all.txt", delim=" ", col_names=F)
names(anc) <- anc_header

# load TOPMed data by population
eur_sites <- read_tsv("/mnt/norbert/data/topmed/singletons/eur/freeze5.singletons.sort.txt2")
afr_sites <- read_tsv("/mnt/norbert/data/topmed/singletons/afr/freeze5.singletons.sort.txt2")

# idsites <- rbind(eur_sites, afr_sites) %>%
#   mutate(BIN=paste0(CHR, ".", ceiling(POS/1e6)))

# idsites <- idsites %>%
idsites <- rbind(eur_sites, afr_sites) %>%
  mutate(BIN=paste0(CHR, ".", ceiling(POS/1e6))) %>%
  group_by(ID) %>%
  mutate(D2n=lead(D2, n = 1L, default=1e6), Dmin=pmin(D2,D2n)) %>%
  dplyr::filter(Dmin>0 & Dmin<5e6)

rm(eur_sites)
rm(afr_sites)
gc()

# count singletons per sample
id_counts <- idsites %>% 
  group_by(ID) %>% 
  summarise(tot=n())

##############################################################################
# Estimate exponential mixture params for inter-mutation distance [IMD]
#
# takes ~5 mins on 8 cores for 10598 samples (~1M sites)
# takes ~2 mins on 8 cores for 1000 samples (~350K sites)
# takes ~2 mins on 8 cores for 1000 samples (~1M sites)
# takes ~3.5 mins on 8 cores for 5000 samples (~5M sites)
# takes ~10 mins on 8 cores for 1000 samples (~2.3M sites)
##############################################################################

#-----------------------------------------------------------------------------
# Run on singleton data
#-----------------------------------------------------------------------------
cluster <- create_cluster(cores = 8)

# deprecated freeze 3 data
# exp_fits_fh <- "/mnt/norbert/data/topmed/fitmix/exp_fits.txt"

# output from 4-component mode
# exp_fits_fh <- "/mnt/norbert/data/topmed/fitmix/exp_fits_freeze5_min.txt"

# output from 4-component minimum distance
exp_fits_fh <- "/mnt/norbert/data/topmed/fitmix/exp_fits_freeze5_min.txt"

# output from running in iteration mode
exp_fits_fh_it <- "/mnt/norbert/data/topmed/fitmix/exp_fits_freeze5_it.txt"
exp_fits_it <- read_tsv(exp_fits_fh_it)

override <- FALSE
if(file.exists(exp_fits_fh) & override==FALSE){
  exp_fits <- read_tsv(exp_fits_fh)  
} else {

  # testing without parallelization
  # exp_fits <- idsites %>%
  #   dplyr::filter(ID %in% testids[1:100]) %>%
  #   group_by(ID) %>%
  #   do(fitExpMix(.$D2, 10)) %>%
  #   as_tibble()
  
  # Initialize cluster and prep data
  # cluster <- create_cluster(cores = 8)
  
  cat("Prepping data...")
  
  run_sites <- merge(idsites, id_counts, by="ID") %>%
    # dplyr::filter(tot>10000) %>%
    # dplyr::filter(D2>0 & D2<2e6) %>%
    # group_by(ID) %>%
    # sample_n(5000) %>%
    ungroup() %>% 
    multidplyr::partition(ID, cluster=cluster)
  
  cat("Complete!")
  
  cat("Prepping cluster...")
  
  run_sites %>%
    cluster_library("tidyverse") %>%
    cluster_library("mixtools") %>%
    cluster_library("MASS") %>%
    cluster_assign_value("getExpLogLik", getExpLogLik) %>%
    cluster_assign_value("getNR2", getNR2) %>%
    cluster_assign_value("fitExpMix", fitExpMix) %>%
    cluster_assign_value("runExpMix", runExpMix)
  
  cat("Complete!")
  
  cat("Running models...")
  
  exp_fits <- run_sites %>%
    # do(fitExpMix(x=.$D2, scale=10, mincomp=3, maxcomp=5, iterate=TRUE)) %>%
    do(fitExpMix(x=.$Dmin, scale=10, mincomp=3, maxcomp=4, iterate=FALSE)) %>%
    as_tibble()
  
  cat("Complete!")
  
  write_tsv(exp_fits, exp_fits_fh, col_names=TRUE)
}

parallel::stopCluster(cluster)

#-----------------------------------------------------------------------------
# get per-individual rates and per-population median rates
#-----------------------------------------------------------------------------
exp_fits2_anc <- exp_fits %>%
  left_join(id_counts, by="ID") %>%
  left_join(anc, by="ID") %>%
  gather(pop, pct_anc, AFR:ME) %>% arrange(ID) %>% #head
  group_by(pop) %>%
  mutate(ntile=ntile(pct_anc, 4)) %>%
  dplyr::filter(pct_anc>0.5) %>%
  dplyr::filter(pop %in% c("AFR", "EUR")) 

exp_fits2_rate_summ <- exp_fits2_anc %>% 
  group_by(param, pop) %>% 
  do(tidy(summary(.$rate)))

exp_fits2_lambda_summ <- exp_fits2_anc %>% 
  group_by(param, pop) %>% 
  do(tidy(summary(.$lambda)))

# get ancestral population for each ID
popids <- exp_fits2_anc %>% 
  dplyr::select(ID, pop) %>% 
  group_by(ID) %>% 
  dplyr::slice(1)

# plot global ancestries of each subsample
left_join(popids, anc, by="ID") %>% 
  gather(anc, prop, AFR:ME) %>%
  arrange(pop, anc, prop) %>%
  ggplot(aes(x=ID, y=prop, fill=anc))+
    geom_bar(stat="identity")+
    facet_wrap(~pop, scales="free_x")+
    scale_fill_brewer(palette="Set1")+
    theme_bw()+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

# subsample popids (n per pop) for testing
subids <- popids %>% group_by(pop) %>% sample_n(100)

# data frame with 4 rate estimates per ind
sample_rates <- exp_fits %>% 
  dplyr::select(ID, param, rate) %>% 
  group_by(ID) %>% 
  spread(param, rate)

# data frame with 4 rate estimates per pop
avg_rates <- exp_fits2_anc %>% 
  group_by(param, pop) %>% 
  summarise(rate=median(rate)) %>% 
  ungroup() %>% 
  mutate(param=paste0(param, "a")) %>% 
  group_by(pop) %>% 
  spread(param, rate)

avg_lambda <- exp_fits2_anc %>% 
  group_by(param, pop) %>% 
  summarise(lambda=median(lambda)) %>% 
  ungroup() %>% 
  mutate(param=paste0(param, "a")) %>% 
  group_by(pop) %>% 
  spread(param, lambda)

# exp_fits2_anc %>% 
#   group_by(param, pop) %>% 
#   summarise(rate=median(rate), lambda=median(lambda)) %>% 
#   ungroup() #%>% 
#   # mutate(param=paste0(param, "a")) %>% 
#   # group_by(pop) %>% 
#   # spread(param, rate)

#-----------------------------------------------------------------------------  
# get per-site component assignments
#-----------------------------------------------------------------------------
assignCluster2 <- function(dist, p1,p2,p3,p4){
  test_rates <- c(p1,p2,p3,p4)
  match <- which.max(lapply(test_rates, function(x) pexp(dist+1, 1/x)-pexp(dist-1,1/x)))
  return(paste0("c", match))
}

test_sites <- idsites %>% 
  left_join(popids, by="ID") %>%
  left_join(sample_rates, by="ID") %>% 
  left_join(avg_rates, by="pop") %>% 
  group_by(ID) %>%
  # dplyr::filter(D2>0) %>%
  # mutate(D2n=lead(D2, n = 1L, default=1e6), Dmin=pmin(D2,D2n)) %>%
  ungroup() %>% 
  mutate(rate_clust=pmap_chr(list(Dmin, p1, p2, p3, p4), assignCluster2)) %>% 
  ungroup() %>%
  dplyr::select(-c(p1:p4, p1a:p4a))

#-----------------------------------------------------------------------------  
# rate vs lambda scatterplots
#-----------------------------------------------------------------------------  
source("ERV_mutation_hotspots/scripts/rate_vs_lambda_scatter.R")

#-----------------------------------------------------------------------------  
# singleton mutation spectra histograms
#-----------------------------------------------------------------------------  
source("ERV_mutation_hotspots/scripts/spectra_plots.R")

#-----------------------------------------------------------------------------  
# compare observed & expected/simulated inter-singleton distance distributions
#-----------------------------------------------------------------------------  
source("ERV_mutation_hotspots/scripts/simulated_distributions.R")

#-----------------------------------------------------------------------------  
# analyze genome-wide distribution of each mixture component 
#-----------------------------------------------------------------------------  
source("ERV_mutation_hotspots/scripts/clusters_by_region.R")

test_sites %>%
  group_by(pop) %>%
  dplyr::summarise(n=n(), count20kb=sum(Dmin<20000), prop=count20kb/n)

#-----------------------------------------------------------------------------  
# analyze DNM data
#-----------------------------------------------------------------------------  
source("ERV_mutation_hotspots/scripts/analyze_clusters_dnms.R")

NWD362844_fits <- idsites %>%
  dplyr::filter(ID=="NWD362844") %>%
  dplyr::filter(D2>0) %>%
  # group_by(ID) %>%
  do(fitExpMix(.$D2, scale=10, mincomp=3, maxcomp=4, iterate=FALSE)) %>%
  as_tibble()

test_sites %>%
group_by(pop, rate_clust) %>%
dplyr::summarise(n=n()) %>%
group_by(pop) %>% mutate(tot=sum(n), prop=n/tot)

