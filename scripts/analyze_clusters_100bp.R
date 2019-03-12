#-----------------------------------------------------------------------------
# NOT RUN: Repeat mixture deconvolution on singletons <=100bp apart
#-----------------------------------------------------------------------------
run_sites <- merge(idsites, id_counts, by="ID") %>%
  # dplyr::filter(tot>10000) %>%
  dplyr::filter(D2>0 & D2<=10) %>%
  # group_by(ID) %>%
  # sample_n(5000) %>%
  ungroup() %>% 
  multidplyr::partition(ID, cluster=cluster)

run_sites %>%
  cluster_library("tidyverse") %>%
  cluster_library("mixtools") %>%
  cluster_library("MASS") %>%
  cluster_assign_value("getExpLogLik", getExpLogLik) %>%
  cluster_assign_value("getNR2", getNR2) %>%
  cluster_assign_value("fitExpMix", fitExpMix) %>%
  cluster_assign_value("runExpMix", runExpMix)

exp_fits_tls <- run_sites %>%
  # do(fitExpMix(x=.$D2, scale=10, mincomp=3, maxcomp=5, iterate=TRUE)) %>%
  do(fitExpMix(x=.$D2, scale=10, mincomp=2, maxcomp=4, iterate=TRUE)) %>%
  as_tibble()

exp_fits_tls %>%
  left_join(popids, by="ID") %>%
  ggplot(aes(x=rate, y=lambda, shape=param, colour=pop))+
  scale_fill_gradient(low="grey50", high="white")+
  # stat_density2d(geom = "raster", aes(fill = ..density..), contour = FALSE)+
  # geom_point(aes(colour=.cluster), size=3, alpha=0.2)+
  geom_point(size=3, alpha=0.6)+
  # facet_wrap(~pop)+
  # facet_grid(pop~ntile)+
  scale_colour_viridis("Ancestry", discrete=TRUE, begin=0, end=0.8)+
  scale_shape_manual("Mixture component", values=c(3, 1, 7, 6))+
  # scale_x_continuous(breaks=c(0:5), labels=10^c(0:5))+
  # scale_y_continuous(breaks=seq(-3, 0, 1), labels=10^seq(-3, 0, 1))+
  scale_x_log10(expand=c(0,0), 
                breaks=c(0,1,10,100),
                labels=c("0", "1", "10", "100"))+
  # scale_y_log10(expand=c(0,0),
  #               limits=c(0.1,1.2),
  #               breaks=c(0.001, 0.01, 0.1, 1))+
  xlab("1/rate")+
  ylab("component contribution (lambda)")+
  annotation_logticks()+
  theme_bw()

sample_rates_tls <- exp_fits_tls %>% 
  dplyr::select(ID, param, rate) %>% 
  group_by(ID) %>% 
  spread(param, rate)


assignCluster2_tls <- function(dist, p1,p2){
  test_rates <- c(p1,p2)
  match <- which.max(lapply(test_rates, function(x) pexp(dist+1, 1/x)-pexp(dist-1,1/x)))
  return(paste0("a", match))
}

test_sites_tls <- idsites %>% 
  dplyr::filter(D2>0 & D2<=100) %>%
  # dplyr::filter(ID %in% subids$ID) %>%
  # left_join(subids, by="ID") %>% #head
  left_join(popids, by="ID") %>%
  left_join(sample_rates_tls, by="ID") %>% #head
  # left_join(avg_rates, by="pop") %>% #head
  # head() %>%
  # dplyr::filter(ID=="NWD100314") %>%
  group_by(ID) %>%
  # dplyr::filter(D2>0) %>%
  mutate(D2n=lead(D2, n = 1L, default=1e6), Dmin=pmin(D2,D2n)) %>% #head
  # mutate(D2n=lead(D2, n = 1L, default=1e6)) %>% #head
  ungroup() %>% #head
  # group_by(CHR, POS) %>%
  mutate(rate_clust=pmap_chr(list(Dmin, p1, p2), assignCluster2_tls)) %>% #,
  # rate_clust_m=map(assignCluster2(D2, D2n, c(p1a, p2a, p3a, p4a)))) %>% #head(20)
  # rowwise() %>%
  # mutate(rate_clust=assignCluster(D2, D2n, c(p1, p2, p3, p4)),
  #        rate_clust_m=assignCluster(D2, D2n, c(p1a, p2a, p3a, p4a))) %>% #head(20)
  # unnest(rate_clust) %>%
  ungroup() %>%
  dplyr::select(-c(p1:p2))

bind_rows(test_sites, test_sites_tls) %>% #dplyr::filter(TYPE=="C_G")
  mutate(motif=substr(MOTIF,3,5)) %>%
  group_by(TYPE, motif, rate_clust, pop) %>%
  # group_by(TYPE, rate_clust, pop) %>%
  summarise(n=n()) %>%
  group_by(rate_clust, pop) %>% #head
  mutate(prop=n/sum(n)) %>%
  ungroup() %>%
  mutate(TYPE=gsub("_", ">", TYPE)) %>%
  ggplot(aes(x=motif, y=prop, fill=rate_clust, alpha=pop))+
  # ggplot(aes(x=TYPE, y=prop, fill=rate_clust))+
  geom_bar(stat="identity", position="dodge")+
  scale_fill_manual(values=mpalette[c(1,2,4,3,5,6)])+
  scale_alpha_manual(values=c(0.3, 1))+
  # facet_wrap(~TYPE, scales="free_x", nrow=1)+
  facet_grid(pop~TYPE, scales="free_x")+
  guides(fill=guide_legend("Mixture component"))+
  ylab("Proportion of singletons")+
  # facet_grid(pop~rate_clust)+
  # facet_wrap(~rate_clust, ncol=1, strip.position="right")+
  theme_bw()+
  theme(strip.text=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.y=element_text(size=12),
        axis.title.x=element_blank(), 
        axis.text.x=element_text(angle=90),
        axis.ticks.x=element_blank(),
        legend.title=element_text(size=16),
        legend.text=element_text(size=12),
        legend.position="bottom")