#   do({ tmp <- with(rle(.$rate_clust=="c2"), lengths[values])
#   # do({ tmp <- with(rle(.$rate_clust=="c1"), sum(lengths > 3))
#       data.frame(Count = tmp)})
# 
# test_sites %>%
#   ungroup() %>%
#   dplyr::filter(ID=="NWD101688") %>%
#   mutate(rate_clust_1=lag(rate_clust, n = 1L, default="c4")) %>% #head
#   mutate(rate_clust_2=lag(rate_clust, n = 2L, default="c4")) %>% #head
#   mutate(rate_clust_3=lag(rate_clust, n = 3L, default="c4")) %>% #head
#   mutate(rate_clust_4=lag(rate_clust, n = 4L, default="c4")) %>% head %>% data.frame
# 
#   mutate(match_1=ifelse(rate_clust==rate_clust_1, TRUE, FALSE)) %>% #head
#   mutate(match_2=ifelse(match_1,
#     ifelse(rate_clust==rate_clust_2, TRUE, FALSE), FALSE)) %>% #head
#   mutate(match_3=ifelse(match_2 & match_1,
#     ifelse(rate_clust==rate_clust_3, TRUE, FALSE), FALSE)) %>% #head
#   mutate(match_4=ifelse(match_3 & match_2 & match_1,
#     ifelse(rate_clust==rate_clust_4, TRUE, FALSE), FALSE)) %>% #head
#   group_by(rate_clust) %>% #head
#   summarise(n=n(), n_match_1=sum(match_1), n_match_2=sum(match_2), n_match_3=sum(match_3), n_match_4=sum(match_4)) %>% 
#   mutate(prop_match=(n_match_1+n_match_2+n_match_3+n_match_4)/n)

##############################################################################
# Calculate autocorrelation between distances
##############################################################################

# get pacf data frame
ac1 <- test_sites %>% 
  ungroup() %>% 
  group_by(ID, pop) %>% 
  # dplyr::select(D2) %>%
  # mutate(D2=log10(D2)) %>% 
  do(tidy(Pacf(.$POS, plot=F)))

# boxplots of acf per pop
ac1 %>% 
  ggplot(aes(x=factor(lag), y=acf, fill=pop))+
  geom_boxplot()+
  # facet_wrap(~rate_clust)+
  xlab("Lag")+
  ylab("Partial autocorrelation")+
  theme_bw()

ac1 %>%
  group_by(lag) %>%
  do(tidy(t.test(.$acf~.$pop)))

# barplot with mean acf per pop
ac1 %>% 
  group_by(lag, pop) %>% 
  summarise(acf=mean(acf)) %>% 
  ggplot(aes(x=lag, y=acf, fill=pop))+
  geom_col(position="dodge")+
  xlab("Lag")+
  ylab("Partial autocorrelation")+
  theme_bw()

ac2 <- test_sites %>% 
  ungroup() %>% 
  group_by(ID) %>% 
  dplyr::select(rate_clust) %>% 
  mutate(rate_clust=as.numeric(gsub("c", "", rate_clust))) %>%
  # mutate(D2=log10(D2)) %>% 
  do(tidy(Pacf(.$rate_clust, plot=F)))

ac2 %>% 
  ggplot(aes(x=lag, y=acf, group=lag))+
  geom_boxplot()+
  xlab("Lag")+
  ylab("Partial autocorrelation")+
  theme_bw()


test_sites %>%
  ungroup() %>%
  dplyr::filter(cl_ID!="UC") %>% 
  dplyr::filter(cl_LEN==cl_NUM) %>% 
  group_by(pop, ID, cl_ID) %>% 
  slice(1L) %>% head
summarise(n=n())

test_sites %>% 
  ungroup() %>% 
  dplyr::filter(cl_ID!="UC") %>% 
  dplyr::filter(cl_LEN==cl_NUM) %>% 
  group_by(pop, TYPE) %>% 
  summarise(n=n()) %>% 
  mutate(tot=sum(n), prop=n/tot) %>% 
  group_by(TYPE) %>% do(tidy(prop.test(.$n,.$tot))) %>% data.frame



##############################################################################
# simulate IMD from per-individual single-parameter exponentials
##############################################################################
exp_fits_wide <- exp_fits %>% 
  as_tibble() %>% 
  gather(variable, value, -(ID:param)) %>% 
  unite(temp, param, variable) %>% 
  spread(temp, value) %>% 
  arrange(ID)

exp_fits_wide <- merge(exp_fits_wide, id_counts, by="ID")

simlens <- unlist(lapply(1:nrow(exp_fits_wide), function(x) simDists(exp_fits_wide, x)))
simlens.df <- data.frame(val=simlens, gp="sim")
simlens.df %>% 
  dplyr::filter(val>=1 & val<2e6) %>%
  ggplot(aes(x=val, fill=gp)) +
  geom_histogram(aes(y = ..density..), position = 'identity', alpha=0.5, bins=1000)+
  xlab("Distance from previous singleton in same individual")+
  theme_bw()+
  theme(axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        strip.text=element_text(size=16))
ggsave("ERV_mutation_hotspots/figs/chr22.D2.sim.hist.png", width=12, height=8)

##############################################################################
# histogram of empirical inter-mutation distance [IMD]
##############################################################################
sites %>% dplyr::filter(D2>=1 & D2<250) %>%
  ggplot(aes(x=D2)) +
  geom_histogram(stat="count", bins=50)+
  scale_y_log10(expand=c(0,0)) +
  xlab("Distance from previous singleton in same individual")+
  theme_bw()+
  theme(axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        strip.text=element_text(size=16),
        legend.position="none")

ggsave("ERV_mutation_hotspots/figs/chr22.D2.hist.png", width=12, height=8)

##############################################################################
# Compare empirical distribution with fitted assuming single parameter
##############################################################################

simexp_e = rexp(1163461,1/352000)
simexp_o = rexp(1163461,1/337200)

cbind(sites, simexp_e, simexp_o) %>% 
  dplyr::select(observed=D2, simulated_e=simexp_e, simulated_o=simexp_o) %>%
  gather(gp, val) %>%
  dplyr::filter(val>=1 & val<2e6) %>%
  # sites %>% dplyr::filter(D1>=2 & D1<100) %>%
  # ggplot(aes(x=factor(D1))) +
  ggplot(aes(x=val, fill=gp)) +
  geom_histogram(aes(y = ..density..), position = 'identity', alpha=0.5, bins=1000)+
  xlab("Distance from previous singleton in same individual")+
  theme_bw()+
  theme(axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        strip.text=element_text(size=16))

ggsave("ERV_mutation_hotspots/figs/chr22.D2.full.hist.sim.png", width=12, height=8)

##############################################################################
# 
##############################################################################
# idsites <- 
sites %>%
  # dplyr::filter(cl_ID!="UC") %>% #dim
  # mutate(cl_CLASS=ifelse(cl_LEN==0, "unclustered",
  #                        ifelse(cl_LEN<10, "1-10",
  #                               ifelse(cl_LEN<100, "11-100", "101+")))) %>% #head
  # mutate(cl_CLASS=factor(cl_CLASS, levels=c("1-10", "11-100", "101+", "unclustered"))) %>%
  mutate(cl_LEN2=ceiling(cl_LEN/5)) %>%
  group_by(cl_LEN2, TYPE) %>%
  summarise(n=n()) %>% #head
  mutate(prop=n/sum(n)) %>% #head
  dplyr::filter(cl_LEN2<50) %>%
  ggplot(aes(x=cl_LEN2, y=prop, fill=TYPE))+
  geom_bar(stat="identity", alpha=0.5, position="identity")+
  facet_wrap(~TYPE)+
  # scale_x_log10()+
  theme_bw()


idsites %>%
  # dplyr::filter(cl_ID!="UC") %>% #dim
  mutate(cl_CLASS=ifelse(cl_LEN==0, "unclustered",
                         ifelse(cl_LEN<10, "1-10",
                                ifelse(cl_LEN<100, "11-100", "101+")))) %>% #head
  mutate(cl_CLASS=factor(cl_CLASS, levels=c("1-10", "11-100", "101+", "unclustered"))) %>%
  group_by(cl_CLASS, TYPE) %>%
  summarise(n=n()) %>% #head
  mutate(prop=n/sum(n)) %>%
  ggplot(aes(x=TYPE, y=prop, fill=TYPE))+
  geom_col()+
  facet_wrap(~cl_CLASS, nrow=1)+
  theme_bw()


##############################################################################
# Get count and exp. param estimate of IMD per 1Mb window
##############################################################################
window_rates <- sites %>%
  dplyr::filter(D2>0) %>%
  mutate(BIN=ceiling(POS/1e5)) %>%
  mutate(cl_len=ifelse(D2==1, "tandem",
                       ifelse(D2<10, "<10bp",
                              ifelse(D2<100, "10-100bp", #">100bp")
                                     ifelse(D2<1000, "100-1000bp", ">1000bp")
                                     #   ifelse(D2<10000, "1000-10000bp", "10000+bp")
                                     # )
                              )
                       )
  )
  ) %>%
  group_by(CHR, BIN, cl_len) %>%
  # group_by(CHR, BIN) %>%
  # dplyr::filter(D2<2e6 & D1<1e6 & D2>10000) %>%
  # dplyr::filter(D2<10000) %>%
  # summarise(fitted_mu_gp=1/getExpParam(D1)*3560, fitted_mu_ind=1/getExpParam(D2)) %>%
  dplyr::summarise(n=n(), fitted_mu_ind=1/getExpParam(D2)) %>%
  # dplyr::filter(n>10000) %>%
  gather("gp", "val", c("n", "fitted_mu_ind")) #%>%
# dplyr::filter(gp=="fitted_mu_ind") #%>%
# group_by(cl_len) %>%
# mutate(val = (val - mean(val)) / sd(val))

total_window_counts <- window_rates %>% 
  dplyr::filter(gp=="n") %>% 
  group_by(CHR, BIN) %>% 
  summarise(n=sum(val)) %>% 
  arrange(n)

cl_levs <- c("tandem", "<10bp", "10-100bp", "100-1000bp", ">1000bp")

sites %>% dplyr::filter(D2>0 & D2<2) %>% summarise(n=n(), prop=n/1151709)

##############################################################################
# plot per-window IMD_s
##############################################################################
merge(window_rates, total_window_counts, by=c("CHR", "BIN")) %>% 
  mutate(cl_len=factor(cl_len, levels=cl_levs)) %>%
  dplyr::filter(n>1500) %>%
  # dplyr::filter(gp=="fitted_mu_ind") %>%
  dplyr::filter(gp=="n") %>%
  group_by(cl_len) %>%
  # dplyr::filter(cl_len %in% c("10-1000bp", "10000+bp")) %>%
  # dplyr::filter(val>10) %>% #head
  # ungroup() %>%
  # group_by(cl_len) %>%
  mutate(z=(val-mean(val))/sd(val)) %>% #head
  # dplyr::filter(abs(z)<8) %>%
  # ggplot(aes(x=BIN, y=val, colour=gp, group=gp))+
  ggplot(aes(x=BIN, y=val, colour=cl_len, group=cl_len))+
  # ggplot(aes(x=BIN, y=z))+
  geom_point(alpha=0.8, size=3)+
  scale_colour_viridis(discrete=T)+
  # scale_colour_manual(values=wes_palette("Darjeeling"))+
  # geom_smooth(span=0.05, se=F)+
  # geom_smooth(span=0.1)+
  # geom_line()+
  # facet_wrap(~CHR, ncol=2, dir="v", strip.position="left")+
  facet_grid(cl_len~CHR, scales="free")+
  # geom_hline(yintercept=0, linetype="dashed")+
  # scale_y_log10()+
  # scale_y_continuous(limits=c(3e-6,1e-5))
  xlab("Window (100kb)")+
  ylab("Count")+
  theme_bw()+
  theme(axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        strip.text=element_text(size=16),
        legend.position="none")

ggsave("ERV_mutation_hotspots/figs/chr22.cluster.counts.png", width=12, height=8)

##############################################################################
# plot histogram of single-component exponential means
##############################################################################
ind_fit <- sites %>% group_by(ID) %>% summarise(exp=1/getExpParam(D2), n=n())

ggplot(ind_fit, aes(x=n))+
  geom_histogram(bins=1000)+
  xlab("# ERVs per individual")+
  theme_bw()+
  theme(axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        strip.text=element_text(size=16),
        legend.position="none")
ggsave("ERV_mutation_hotspots/figs/chr22.ind.count.hist.png", width=12, height=8)


window_rates %>% 
  dplyr::filter(val<300000) %>%
  spread(gp, val) %>%
  mutate(v2 = fitted_mu_gp/fitted_mu_ind) %>%
  dplyr::filter(v2<1.2 & v2>0.8) %>%
  ggplot(aes(x=BIN, y=v2))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~CHR, ncol=2)+
  # scale_y_continuous(limits=c(3e-6,1e-5))
  theme_bw()

window_rates %>% 
  dplyr::filter(val<300000 & val>150000) %>%
  ungroup() %>%
  spread(gp, val) %>% na.omit() %>% #head
  summarise(cor=cor(fitted_mu_gp, fitted_mu_ind, method="spearman"))

pred1mb <- read_tsv("/mnt/norbert/data/bridges/predicted_1Mb.txt")
names(pred1mb) <- c("file", "CHR", "BIN", "C1", "C2")

pred1mb_notype <- pred1mb %>% 
  group_by(CHR, BIN) %>% 
  summarise(C1=sum(C1), C2=sum(C2)) %>% 
  mutate(exp=C1+C2)

window_counts <- sites %>% group_by(CHR, BIN) %>% summarise(obs=n())
window_counts_merged <- merge(window_counts, pred1mb_notype, by=c("CHR", "BIN")) %>% 
  mutate(resid=obs-exp)

mean_dist <- window_rates %>% 
  dplyr::filter(gp=="n") %>%
  # dplyr::filter(gp=="fitted_mu_ind") %>%
  dplyr::select(-gp) %>% 
  spread(cl_len, val) 

window_counts_merged2 <- merge(window_counts_merged, mean_dist, by=c("CHR", "BIN")) %>% 
  dplyr::filter(obs>8000)


m1 <- lm(resid~`<10bp`+`10-1000bp`+`1000-10000bp`+`10000+bp`, data=window_counts_merged2)

summary(m1)

af <- anova(m1)
afss <- af$"Sum Sq"
print(cbind(af,PctExp=afss/sum(afss)*100))


#-----------------------------------------------------------------------------
# initial try at correlation between ext branch length and % clustered
# assuming ext branches are exponentially distributed (not correct)
# --revised simulation with proper demographic model in separate jupyter nb
#-----------------------------------------------------------------------------

n_branches <- 10000
chr_size <- 1e7 # 10Mb
mu <- 1.5e-8
mean_gens <- 200
ext_branch_lengths <- round(rexp(n_branches,1/mean_gens))

full_data <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(full_data) <- c("branch_length", "interval")

for(i in 1:n_branches){
  branch_length <- ext_branch_lengths[i]
  n_mutations <- rpois(1, branch_length*chr_size*mu)
  sites <- sort(round(runif(n_mutations, 1, chr_size)))
  intervals <- sites[-1]-sites[-n_mutations]
  n_intervals <- length(intervals)
  branch_data <- data.frame(branch_length=rep(branch_length, n_intervals), interval=intervals)
  
  full_data <- bind_rows(full_data, branch_data)
  
}

plot_data <- full_data %>% 
  mutate(cl=ifelse(interval<20000, TRUE, FALSE)) %>% #summarise(n=n(), ncl=sum(cl), prop=ncl/n)
  group_by(branch_length) %>% 
  summarise(n=n(), ncl=sum(cl)) %>% #head(20)
  mutate(propcl=ncl/n)

# plotdata$branch_length, plotdata$propcl

plot_data %>% #summarise(meanprop=mean(propcl))
  ggplot(aes(x=branch_length, y=propcl))+
  geom_point(alpha=0.5)+
  geom_smooth(method="lm")+
  geom_vline(xintercept=mean_gens, linetype="dashed")+
  scale_x_log10()+
  scale_y_log10()+
  annotation_logticks()+
  xlab("External branch length (generations)")+
  ylab("Proportion of intervals <20kb apart")+
  ggtitle("[mean branch length: 200; mu=1.5e-8]")+
  theme_bw()
