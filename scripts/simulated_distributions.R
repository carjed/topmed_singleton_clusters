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
ggsave(paste0(scriptdir, "/figs/distance_distributions.png"), width=8, height=4)


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
sim_sites <- read_tsv(paste0(scriptdir, "/data/sim_dists_dip.txt"))
sim_sites_rc <- read_tsv(paste0(scriptdir, "/data/sim_dists_dip_rc.txt"))

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
