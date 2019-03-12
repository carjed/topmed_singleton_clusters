#-----------------------------------------------------------------------------
# BRIDGES data
#-----------------------------------------------------------------------------
br_sites <- read_tsv("/mnt/norbert/data/bridges/bridges_ervs_dist_sort_5M.txt", col_names=F)
header <- "CHR,POS,ALT,REF,AA,AN,MOTIF,TYPE,Category2,BIN,MASK,S,ID,D1,D2"
names(br_sites) <- as.list(strsplit(header,","))[[1]]

# br_testids <- head(unique(br_sites$ID), 100)

br_sites <- br_sites %>%
  dplyr::select(CHR, POS, REF, ALT, TYPE, MOTIF, ID, D2)

#-----------------------------------------------------------------------------
# Get mixture parameters for BRIDGES data
#-----------------------------------------------------------------------------

br_id_counts <- br_sites %>% 
  group_by(ID) %>% 
  summarise(tot=n())

run_sites <- merge(br_sites, br_id_counts, by="ID") %>%
  dplyr::filter(tot>5000 & tot<12000) %>%
  group_by(ID) %>%
  do(sample_n(.,5000)) %>%
  dplyr::filter(D2>0 & D2<2e6) %>%
  ungroup()

br_sites_part <- multidplyr::partition(run_sites, ID, cluster=cluster)
br_sites_part %>%
  cluster_library("tidyverse") %>%
  cluster_library("mixtools") %>%
  cluster_library("MASS") %>%
  cluster_assign_value("getExpLogLik", getExpLogLik) %>%
  cluster_assign_value("getNR2", getNR2) %>%
  cluster_assign_value("fitExpMix", fitExpMix) %>%
  cluster_assign_value("runExpMix", runExpMix)

exp_fits_br <- br_sites_part %>%
  do(fitExpMix(x=.$D2, scale=10, mincomp=3, maxcomp=4)) %>%
  as_tibble()

exp_fits_br2 <- exp_fits_br %>%
  mutate(tot=5000, pop="EUR (BRIDGES)", pct_anc=1, ntile=4)


exp_fits2b <- merge(exp_fits_b, bin_counts, by="BIN")

# rbind(exp_fits2_anc, exp_fits_br2) %>%  
#   as_tibble() %>% #head
#   dplyr::filter(lambda>0.0001) %>%
#   ggplot(aes(x=rate, y=lambda))+
#     scale_fill_gradient(low="grey50", high="white")+
#     stat_density2d(geom = "raster", aes(fill = ..density..), contour = FALSE)+
#     # geom_point(aes(colour=.cluster), size=3, alpha=0.2)+
#     geom_point(aes(colour=param), size=3, alpha=0.2)+
#     facet_wrap(~pop)+
#     # facet_grid(pop~ntile)+
#     scale_colour_viridis(discrete=TRUE)+
#     scale_x_log10(expand=c(0,0), 
#                   breaks=c(0,1,10,100,1000,10000,100000,1e6),
#                   labels=c("0", "1", "10", "100", "1,000", "10,000", "100,000", "1M"))+
#     scale_y_log10(expand=c(0,0),
#                   limits=c(0.001,1.2),
#                   breaks=c(0.01, 0.1, 1))+
#     annotation_logticks()+
#     theme_bw()