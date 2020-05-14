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

# load yaml package to parse args
install.packages("yaml", quiet=TRUE)
library(yaml)

# parse args
args <- yaml.load_file("../rconfig.yaml")
attach(args)

source(paste0(scriptdir, "/cluster_functions.R"))

assignCluster2 <- function(dist, p1,p2,p3,p4){
  test_rates <- c(p1,p2,p3,p4)
  match <- which.max(lapply(test_rates, function(x) pexp(dist+1, 1/x)-pexp(dist-1,1/x)))
  return(paste0("c", match))
}

##############################################################################
# Read and prep data
##############################################################################

id_counts <- tibble(ID=character(), tot=integer())
test_sites <- tibble()

# African, CSAsian, Easian, European, Native American, Oceania, Middle Eastern
anc_header <- c("ID", "AFR", "CSA", "EAS", "EUR", "NAT", "OC", "ME")

# freeze 5 ancestry
anc <- read_delim(paste0(datadir, "/ancestry/global_list_b38all.txt"), delim=" ", col_names=F)
names(anc) <- anc_header

# load TOPMed data by population
cluster <- create_cluster(cores = 10)

pops <- c("afr", "eas", "eur")
for(pop in pops){
  sites <- read_tsv(paste0(datadir, "/singletons/", pop, "/freeze5.singletons.sort.txt2"))
  # eur_sites <- read_tsv(paste0(datadir, "/singletons/eur/freeze5.singletons.sort.txt2"))
  # afr_sites <- read_tsv(paste0(datadir, "/singletons/afr/freeze5.singletons.sort.txt2"))
  # eas_sites <- read_tsv(paste0(datadir, "/singletons/eas/freeze5.singletons.sort.txt2"))

  # idsites <- rbind(eur_sites, afr_sites) %>%
  #   mutate(BIN=paste0(CHR, ".", ceiling(POS/1e6)))

  # idsites <- idsites %>%
  # idsites <- rbind(eur_sites, afr_sites, eas_sites) %>%
  idsites <- sites %>%
    # idsites <- eur_sites %>%
    mutate(BIN=paste0(CHR, ".", ceiling(POS/1e6))) %>%
    group_by(ID) %>%
    mutate(D2n=lead(D2, n = 1L, default=1e6), Dmin=pmin(D2,D2n)) %>%
    dplyr::filter(Dmin>0 & Dmin<5e6)

  rm(sites)
  # rm(eur_sites)
  # rm(afr_sites)
  # rm(eas_sites)
  gc()

  #-----------------------------------------------------------------------------
  # Run on singleton data
  #-----------------------------------------------------------------------------

  # cached output from 4-component minimum distance
  exp_fits_fh <- paste0(datadir, "/fitmix/exp_fits_freeze5_min_", pop, ".txt")

  # output from running in iteration mode
  # exp_fits_fh_it <- paste0(datadir, "/fitmix/exp_fits_freeze5_it.txt")
  # exp_fits_it <- read_tsv(exp_fits_fh_it)

  override <- FALSE
  if(file.exists(exp_fits_fh) & override==FALSE){
    # exp_fits <- read_tsv(exp_fits_fh)

    exp_fits <- file.info(list.files(paste0(datadir, "/fitmix/"), pattern="txt", full.names=TRUE)) %>%
      rownames %>%
      map_dfr(read_tsv)

    if(!exists("exp_fits2_anc")){
      #-----------------------------------------------------------------------------
      # get per-individual rates and per-population median rates
      #-----------------------------------------------------------------------------
      exp_fits2_anc <- exp_fits %>%
        left_join(id_counts, by="ID") %>%
        left_join(anc, by="ID") %>%
        gather(pop, pct_anc, AFR:ME) %>% arrange(ID) %>%
        group_by(pop) %>%
        mutate(ntile=ntile(pct_anc, 4)) %>%
        dplyr::filter(pct_anc>0.5) %>%
        dplyr::filter(pop %in% c("AFR", "EUR", "EAS"))

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
    }

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



  #-----------------------------------------------------------------------------
  # get per-site component assignments
  #-----------------------------------------------------------------------------

  test_sites_sub <- idsites %>%
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

  # count singletons per sample
  id_counts_sub <- idsites %>%
    group_by(ID) %>%
    summarise(tot=n())

  id_counts <- bind_rows(id_counts, id_counts_sub)
  test_sites <- bind_rows(test_sites, test_sites_sub)
}


##############################################################################
# Estimate exponential mixture params for inter-mutation distance [IMD]
#
# takes ~5 mins on 8 cores for 10598 samples (~1M sites)
# takes ~2 mins on 8 cores for 1000 samples (~350K sites)
# takes ~2 mins on 8 cores for 1000 samples (~1M sites)
# takes ~3.5 mins on 8 cores for 5000 samples (~5M sites)
# takes ~10 mins on 8 cores for 1000 samples (~2.3M sites)
##############################################################################



parallel::stopCluster(cluster)


#-----------------------------------------------------------------------------
# get per-individual rates and per-population median rates
#-----------------------------------------------------------------------------
exp_fits2_anc <- exp_fits %>%
  left_join(id_counts, by="ID") %>%
  left_join(anc, by="ID") %>%
  gather(pop, pct_anc, AFR:ME) %>% arrange(ID) %>%
  group_by(pop) %>%
  mutate(ntile=ntile(pct_anc, 4)) %>%
  dplyr::filter(pct_anc>0.5) %>%
  dplyr::filter(pop %in% c("AFR", "EUR", "EAS"))

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

tsc <- test_sites %>% mutate(Dmin=ifelse(Dmin>10000, 10001, Dmin)) %>% group_by(pop, Dmin) %>% count()


tsc %>%
  group_by(pop) %>%
  mutate(tot=sum(n)) %>%
  mutate(prop=n/tot) %>%
  dplyr::filter(Dmin<=100) %>%
  ggplot(aes(x=Dmin, y=prop, colour=pop, group=pop))+
    geom_point(alpha=0.5)+
  geom_line()+
  scale_x_log10()+
  scale_y_log10()+
  scale_colour_manual("Ancestry", values=c("#FF7F00", "#E31A1C", "#33A02C"))+
  xlab("Distance to nearest singleton in same individual (bp)")+
  ylab("Proportion of all singletons")+
  annotation_logticks()+
  theme_bw()+
  theme(#strip.text=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=16),
        axis.text.x=element_text(size=12),
        # axis.ticks.x=element_blank(),
        legend.title=element_text(size=16),
        legend.text=element_text(size=12))

ggsave(paste0(projdir, "/figs/isd_prop.png"), width=12, height=8)

#-----------------------------------------------------------------------------
# rate vs lambda scatterplots
#-----------------------------------------------------------------------------
source(paste0(scriptdir, "/rate_vs_lambda_scatter.R"))

#-----------------------------------------------------------------------------
# singleton mutation spectra histograms
#-----------------------------------------------------------------------------
source(paste0(scriptdir, "/spectra_plots.R"))

#-----------------------------------------------------------------------------
# compare observed & expected/simulated inter-singleton distance distributions
#-----------------------------------------------------------------------------
source(paste0(scriptdir, "/simulated_distributions.R"))

#-----------------------------------------------------------------------------
# analyze genome-wide distribution of each mixture component
#-----------------------------------------------------------------------------
source(paste0(scriptdir, "/clusters_by_region.R"))

test_sites %>%
  group_by(pop) %>%
  dplyr::summarise(n=n(), count20kb=sum(Dmin<20000), prop=count20kb/n)

#-----------------------------------------------------------------------------
# analyze DNM data
#-----------------------------------------------------------------------------
# source(paste0(scriptdir, "/analyze_clusters_dnms.R"))


