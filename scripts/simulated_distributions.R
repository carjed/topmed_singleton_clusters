#-----------------------------------------------------------------------------
#plot distance distributions
#-----------------------------------------------------------------------------

simdist_a <- function(x, y){
  rexp(n=x*100000, rate=y)
}

# DEPRECATED simulated marginal based on overall median
# sim_marginal <- exp_fits2_anc %>% 
#   group_by(pop) %>% 
#   summarise(lambda=median(tot)/100000, rate=median(tot)/3e9) %>%
#   mutate(D2=map2(lambda, rate, simdist_a)) %>% 
#   unnest() %>%
#   mutate(group="exp_marginal")

# simulation based on mixture of 4 fitted components (medians)
sim_mixture <- exp_fits2_anc %>% 
  group_by(pop, param) %>% 
  summarise(lambda=median(lambda), rate=1/median(rate)) %>%
  mutate(D2=map2(lambda, rate, simdist_a)) %>% 
  unnest() %>%
  mutate(group="exp_mixture") %>%
  dplyr::select(pop, D2, group)

# simulation based on mixture of 4k fitted components (4 per sample)
# sim_mixture <- exp_fits2_anc %>% 
#   group_by(pop, param) %>% 
#   summarise(lambda=median(lambda), rate=1/median(rate)) %>%
#   mutate(D2=map2(lambda, rate, simdist_a)) %>% 
#   unnest() %>%
#   mutate(group="exp_mixture") %>%
#   dplyr::select(pop, D2, group)

# EMPIRICAL SIMULATION 1a:
# simulation based on mixture of 1k components (1 per sample)
id_anc <- exp_fits2_anc %>% 
  dplyr::filter(param=="p1") %>% 
  dplyr::select(ID, pop) 

sim_marginal2 <- left_join(id_counts, id_anc, by="ID") %>% 
  group_by(pop) %>% 
  mutate(D2=map2(tot/100000, tot/3e9, simdist_a)) %>% 
  unnest() %>%
  mutate(group="exp_marginal_1a")

test_sites %>% 
  group_by(pop) %>%
  dplyr::filter(D2<3e6) %>%
  sample_n(1e5) %>%
  mutate(group=".observed") %>%
  bind_cols(data.frame(panel=rep(c("4-component", "1-component"), 1e5))) %>% #head
  dplyr::select(pop, Dmin, group, panel) %>%
  bind_rows(sim_mixture %>% mutate(Dmin=D2, panel="4-component")) %>% #head()
  bind_rows(sim_marginal2 %>% sample_n(1e5) %>% mutate(Dmin=D2, panel="1-component")) %>%
  mutate(D2=log(Dmin)) %>%
  ggplot(.)+
  geom_line(aes(x=D2, colour=pop, linetype=group, alpha=group), 
            size=1, stat="density")+
  scale_x_continuous(limits=c(-1,15), 
                     breaks=log(c(1,150,20000,300000)), 
                     labels=c(1,150,20000,300000))+
  scale_colour_manual("Ancestry", values=c("#FF7F00", "#33A02C"))+
  scale_alpha_manual(values=c(1,0.7,0.7))+
  # scale_linetype("Data")+
  # scale_linetype_manual(guide = 'none', values=c("solid", "twodash", "dotted"))+
  facet_wrap(~panel)+
  xlab("Inter-singleton distance")+
  ylab("density")+
  theme_bw()+
  theme(legend.position="bottom")+
  guides(linetype = FALSE, alpha=FALSE)
ggsave("/mnt/norbert/home/jedidiah/projects/ERV_mutation_hotspots/figs/distance_distributions.png", width=8, height=4)



# plot inter-singleton distance distribution of empirical simulation
# test_sites %>% 
#   group_by(pop) %>%
#   sample_n(1e5) %>%
#   dplyr::filter(D2<3e6) %>%
#   mutate(group=".observed") %>%
#   dplyr::select(pop, D2, group) %>%
#   bind_rows(sim_marginal2 %>% sample_n(1e5)) %>%
#   mutate(D2=log(D2)) %>%
#   ggplot(.)+
#   geom_line(aes(x=D2, colour=pop, linetype=group), 
#             size=1, stat="density")+
#   scale_x_continuous(limits=c(-1,15), 
#                      breaks=log(c(1,150,20000,300000)), 
#                      labels=c(1,150,20000,300000))+
#   scale_colour_manual("Ancestry", values=c("#FF7F00", "#33A02C"))+
#   scale_linetype("Data")+
#   xlab("Inter-singleton distance")+
#   ylab("density")+
#   theme_bw()+
#   theme(legend.position="bottom")
# ggsave("/mnt/norbert/home/jedidiah/projects/ERV_mutation_hotspots/figs/distance_distributions2.png", width=8, height=4)

# EMPIRICAL SIMULATION 1b:
# based on mixture of 1k components (1 per sample) with variable mutation rate
simdist_mix <- function(x){
  c(rexp(n=0.9*x, rate=x/3e9), rexp(n=0.1*x, rate=2*x/3e9))
}

sim_marginal2_var_mu <- left_join(id_counts, id_anc, by="ID") %>% 
  group_by(pop) %>%
  mutate(D2=map(tot, simdist_mix)) %>% 
  unnest() %>%
  mutate(group="exp_marginal_1b")

# COALESCENT SIMULATIONS
# load data output from msprime Jupyter notebook
sim_sites <- read_tsv("/mnt/norbert/home/jedidiah/projects/ERV_mutation_hotspots/data/sim_dists_dip.txt")
sim_sites_rc <- read_tsv("/mnt/norbert/home/jedidiah/projects/ERV_mutation_hotspots/data/sim_dists_dip_rc.txt")

# calculate proportion of inter-singleton distances <20kb attributable to 
# external branch length heterogeneity in empirical simulation

# empirical sim, constant mu
sim_marginal2 %>% 
  group_by(pop) %>%
  dplyr::summarise(n=n(), count20kb=sum(D2<8700), prop=count20kb/n)

sim_marginal2 %>% 
  group_by(pop) %>%
  dplyr::summarise(n=n(), count20kb=sum(D2<120000), prop=count20kb/n)

# empirical sim, variable mu (10% of genome subject to 2x higher mu)
sim_marginal2_var_mu %>% 
  group_by(pop) %>%
  dplyr::summarise(n=n(), count20kb=sum(D2<20000), prop=count20kb/n)

# coalescent sim (EUR only!)
sim_sites %>%
  dplyr::summarise(n=n(), count20kb=sum(dist<20000), prop=count20kb/n)

# 
# #-----------------------------------------------------------------------------
# # Run mixture deconvolution on simulated data (coalescent)
# #-----------------------------------------------------------------------------
# id_counts_sim <- sim_sites %>% 
#   group_by(ind_dip) %>% 
#   summarise(tot=n())
# 
# cluster <- create_cluster(cores = 8)
# 
# run_sites_sim <- merge(sim_sites, id_counts_sim, by="ind_dip") %>% 
#   # dplyr::filter(tot>5000) %>%
#   dplyr::filter(dist<2e6) %>%
#   multidplyr::partition(ind_dip, cluster=cluster)
# 
# run_sites_sim %>%
#   cluster_library("tidyverse") %>%
#   cluster_library("mixtools") %>%
#   cluster_library("MASS") %>%
#   cluster_assign_value("getExpLogLik", getExpLogLik) %>%
#   cluster_assign_value("getNR2", getNR2) %>%
#   cluster_assign_value("fitExpMix", fitExpMix) %>%
#   cluster_assign_value("runExpMix", runExpMix)
# 
# exp_fits_sim <- run_sites_sim %>%
#   do(fitExpMix(x=.$dist, scale=10, mincomp=2, maxcomp=2, iterate=TRUE)) %>%
#   as_tibble()
# 
# parallel::stopCluster(cluster)
# 
# write_tsv(exp_fits_sim, "/mnt/norbert/data/topmed/fitmix/exp_fits_sim.txt", col_names=TRUE)
# 
# #-----------------------------------------------------------------------------
# # Plot parameter estimates from simulated data
# #-----------------------------------------------------------------------------
# exp_fits_sim %>% 
#   dplyr::filter(lambda>1e-5) %>% 
#   ggplot(aes(x=rate, y=lambda, shape=param))+
#   scale_fill_gradient(low="grey50", high="white")+
#   # stat_density2d(geom = "raster", aes(fill = ..density..), contour = FALSE)+
#   # geom_point(aes(colour=.cluster), size=3, alpha=0.2)+
#   geom_point(size=3, alpha=0.6)+
#   # facet_wrap(~pop)+
#   # facet_grid(pop~ntile)+
#   # scale_colour_viridis("Ancestry", discrete=TRUE, begin=0, end=0.8)+
#   scale_shape_manual("Mixture component", values=c(3, 1, 7, 6))+
#   # scale_x_continuous(breaks=c(0:5), labels=10^c(0:5))+
#   # scale_y_continuous(breaks=seq(-3, 0, 1), labels=10^seq(-3, 0, 1))+
#   scale_x_log10(expand=c(0,0), 
#                 breaks=c(0,1,10,100,1000,10000,100000,1e6),
#                 labels=c("0", "1", "10", "100", "1,000", "10,000", "100,000", "1M"))+
#   scale_y_log10(expand=c(0,0),
#                 limits=c(0.0001,1.2),
#                 breaks=c(0.001, 0.01, 0.1, 1))+
#   xlab("1/rate")+
#   ylab("component contribution (lambda)")+
#   annotation_logticks()+
#   theme_bw()
# 
# 
# simdist_core <- function(lambda, rate){
#   rexp(lambda*20000, rate)
# }
# 
# core_dists <- avg_rates %>%
#   gather(param, rate, p1a:p4a) %>%
#   left_join(avg_lambda %>% gather(param, lambda, p1a:p4a), by=c("pop", "param")) %>% #head
#   group_by(pop, param) %>%
#   mutate(D2=map2(lambda, 1/rate, simdist_core)) %>% 
#   unnest() %>%
#   mutate(rank=ntile(D2, 100)) %>%
#   dplyr::filter(rank>25 & rank<=75) %>%
#   summarise(max=max(D2), min=min(D2)) %>%
#   gather(limit, value, max:min) %>%
#   mutate(param=gsub("a", "", param)) %>%
#   mutate(key=paste(pop, param, limit, sep="_")) %>%
#   mutate(value=ceiling(value)) %>%
#   ungroup() %>%
#   dplyr::select(key, value) %>%
#   spread(key, value)
# 
# class_counts <- test_sites %>%
#   group_by(pop, rate_clust) %>%
#   count() %>%
#   group_by(pop) %>%
#   mutate(tot=sum(n), prop=n/tot)
# 
# # cluster 1 AFR/EUR
# test_sites %>%
#   group_by(pop) %>%
#   dplyr::summarise(n=n(), 
#                    count20kb=sum(Dmin<core_dists$AFR_p1_max), 
#                    prop=count20kb/n)
# 
# sim_marginal2 %>%
#   group_by(pop) %>%
#   dplyr::summarise(n=n(), 
#                    count20kb=sum(D2<core_dists$AFR_p1_max), 
#                    prop=count20kb/n)
# 
# # cluster 2 AFR
# test_sites %>%
#   dplyr::filter(pop=="AFR") %>%
#   dplyr::summarise(n=n(), 
#                    count20kb=sum(Dmin<core_dists$AFR_p2_max & Dmin>core_dists$AFR_p2_min), 
#                    prop=count20kb/n)
# 
# sim_marginal2 %>%
#   dplyr::filter(pop=="AFR") %>%
#   mutate(Dmin=D2) %>%
#   dplyr::summarise(n=n(), 
#                    count20kb=sum(Dmin<core_dists$AFR_p2_max & Dmin>core_dists$AFR_p2_min), 
#                    prop=count20kb/n)
# 
# # cluster 3 AFR
# test_sites %>%
#   dplyr::filter(pop=="AFR") %>%
#   dplyr::summarise(n=n(), 
#                    count20kb=sum(Dmin<core_dists$AFR_p3_max & Dmin>core_dists$AFR_p3_min), 
#                    prop=count20kb/n)
# 
# sim_marginal2 %>%
#   dplyr::filter(pop=="AFR") %>%
#   mutate(Dmin=D2) %>%
#   dplyr::summarise(n=n(), 
#                    count20kb=sum(Dmin<core_dists$AFR_p3_max & Dmin>core_dists$AFR_p3_min), 
#                    prop=count20kb/n)
# 
# # cluster 2 EUR
# test_sites %>%
#   dplyr::filter(pop=="EUR") %>%
#   dplyr::summarise(n=n(), 
#                    count20kb=sum(Dmin<core_dists$EUR_p2_max & Dmin>core_dists$EUR_p2_min),
#                    prop=count20kb/n)
# 
# sim_marginal2 %>%
#   dplyr::filter(pop=="EUR") %>%
#   mutate(Dmin=D2) %>%
#   dplyr::summarise(n=n(), 
#                    count20kb=sum(Dmin<core_dists$EUR_p2_max & Dmin>core_dists$EUR_p2_min), 
#                    prop=count20kb/n)
# 
# # cluster 3 EUR
# test_sites %>%
#   dplyr::filter(pop=="EUR") %>%
#   dplyr::summarise(n=n(), 
#                    count20kb=sum(Dmin<core_dists$EUR_p3_max & Dmin>core_dists$EUR_p3_min), 
#                    prop=count20kb/n)
# 
# sim_marginal2 %>%
#   dplyr::filter(pop=="EUR") %>%
#   mutate(Dmin=D2) %>%
#   dplyr::summarise(n=n(), 
#                    count20kb=sum(Dmin<core_dists$EUR_p3_max & Dmin>core_dists$EUR_p3_min), 
#                    prop=count20kb/n)
