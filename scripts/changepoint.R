#-----------------------------------------------------------------------------
# TESTING changepoint analysis
#-----------------------------------------------------------------------------

# do(tidy(cpt.np(.$prop, method="PELT", penalty="SIC", minseglen=2, class=FALSE))) %>% data.frame




# idsites %>% 
#   dplyr::filter(cl_ID!="UC") %>% 
#   # dplyr::filter(rate_clust %in% c("c1", "c2")) %>% 
#   mutate(gp=ntile(cl_LEN, 50)) %>% 
#   # mutate(cl_LEN=ifelse(cl_LEN>100 & cl_LEN<1000, ceiling(cl_LEN/10)*10, cl_LEN)) %>%
#   group_by(cl_LEN) %>%
#   group_by(cl_LEN, TYPE) %>% 
#   summarise(n=n()) %>% 
#   # summarise(mean.mpg = mean(mpg, na.rm = TRUE),
#   #           sd.mpg = sd(mpg, na.rm = TRUE),
#   #           n.mpg = n()) %>%
#   # mutate(se.mpg = sd.mpg / sqrt(n.mpg),
#   #        lower.ci.mpg = mean.mpg - qt(1 - (0.05 / 2), n.mpg - 1) * se.mpg,
#   #        upper.ci.mpg = mean.mpg + qt(1 - (0.05 / 2), n.mpg - 1) * se.mpg)
#   group_by(cl_LEN) %>% 
#   mutate(prop=n/sum(n), tot=sum(n), error = sqrt((prop * (1-prop))/n)) %>%
#   dplyr::filter(cl_LEN<=100) %>%
#   ggplot(aes(x=cl_LEN, y=prop, fill=TYPE, group=TYPE))+
#   # scale_x_log10(breaks=sort(c(5*10^(0:3), 10^(0:3))))+
#   # scale_x_log10(expand=c(0,0.02), breaks=sort(c(1:9, seq(10,100,by=10))))+
#   # scale_y_continuous(breaks=seq(0,0.4,by=0.05))+
#   # geom_point(aes(size=tot), alpha=0.5)+
#   # geom_point(alpha=0.5)+
#   geom_bar(stat="identity", alpha=0.5)+
#   geom_errorbar(aes(ymin = prop - 1.96*error, ymax = prop + 1.96*error, color=TYPE),
#                 position = position_dodge(0.9), width=0.5)+
#   facet_wrap(~TYPE, ncol=1, scales="free_y", strip.position="right")+
#   # geom_smooth(span=0.1, se=F)+
#   # geom_line(alpha=0.5)+
#   # annotation_logticks(sides="b")+
#   xlab("distance to nearest singleton")+
#   ylab("singleton fraction")+
#   theme_classic()+
#   theme(legend.position="none")



# 
exp_fits_c1 <- test_sites %>%
  dplyr::filter(rate_clust=="c2") %>%
  dplyr::filter(D2<20000) %>%
  group_by(ID) %>%
  do(fitExpMix(x=.$D2, scale=10, mincomp=2, maxcomp=3, iterate=TRUE)) %>%
  # do(fitExpMix(x=.$D2, scale=10, mincomp=2, maxcomp=2)) %>%
  as_tibble()

#-----------------------------------------------------------------------------
# get frequencies of all possible component sequences by individual variant
#-----------------------------------------------------------------------------
# seqs <- permutations(x=paste("c", 1:4, sep=""), k=5, replace=T) %>% 
seqs <- permutations(x=1:4, k=2, replace=T) %>% 
  data.frame %>% 
  unite(seq, X1:X2, sep="") %>%
  dplyr::select(seq) %>% 
  as.list %>% unlist()

test_seq_counts <- test_sites %>% 
  ungroup() %>% 
  group_by(ID) %>% 
  mutate(rate_clust=gsub("c", "", rate_clust)) %>%
  summarise(seq=paste(rate_clust, collapse="")) %>% 
  group_by(ID) %>% 
  do(data.frame(ID=.$ID, subseq=seqs, x=str_count(.$seq, seqs))) %>%
  group_by(subseq) %>%
  summarise(n=sum(x)) %>%
  arrange(desc(n))

#-----------------------------------------------------------------------------
# get cluster run frequencies
#-----------------------------------------------------------------------------
t1 <- test_sites %>%
  ungroup() %>%
  group_by(ID, pop) %>%
  do(data.frame(unclass(rle(.$rate_clust)))) 

t1 <- idsites %>% 
  dplyr::filter(cl_ID!="UC") %>% 
  # dplyr::filter(rate_clust %in% c("c1", "c2")) %>% 
  mutate(gp=ntile(cl_LEN, 50)) %>% 
  # mutate(cl_LEN=ifelse(cl_LEN>100 & cl_LEN<1000, ceiling(cl_LEN/10)*10, cl_LEN)) %>%
  group_by(cl_LEN) %>%
  group_by(cl_LEN, TYPE) %>% 
  summarise(n=n()) %>% 
  # summarise(mean.mpg = mean(mpg, na.rm = TRUE),
  #           sd.mpg = sd(mpg, na.rm = TRUE),
  #           n.mpg = n()) %>%
  # mutate(se.mpg = sd.mpg / sqrt(n.mpg),
  #        lower.ci.mpg = mean.mpg - qt(1 - (0.05 / 2), n.mpg - 1) * se.mpg,
  #        upper.ci.mpg = mean.mpg + qt(1 - (0.05 / 2), n.mpg - 1) * se.mpg)
  group_by(cl_LEN) %>% 
  mutate(prop=n/sum(n), tot=sum(n), error = sqrt((prop * (1-prop))/n)) %>%
  dplyr::filter(cl_LEN<=100) %>%
  ungroup() %>%
  group_by(TYPE) %>%
  arrange(cl_LEN) %>% #head
  do(tidy(bcp(.$prop)$posterior.prob)) %>%
  mutate(cl_LEN=row_number(), breakpt="BREAK") %>%
  dplyr::filter(x>0.5)

#-----------------------------------------------------------------------------
# get frequencies of all possible component sequences by cluster run
#-----------------------------------------------------------------------------
rc_counts <- test_sites %>% 
  group_by(pop, rate_clust) %>% 
  summarise(n=n()) %>% 
  mutate(prop=n/sum(n)) %>% 
  dplyr::select(-n) %>% #head
  spread(rate_clust, prop) %>%
  mutate(c2_has1=c2/(1-c1), c3_has1=c3/(1-c1), c4_has1=c4/(1-c1),
         c1_has2=c1/(1-c2), c3_has2=c3/(1-c2), c4_has2=c4/(1-c2),
         c1_has3=c1/(1-c3), c2_has3=c2/(1-c3), c4_has3=c4/(1-c3),
         c1_has4=c1/(1-c4), c2_has4=c2/(1-c4), c3_has4=c3/(1-c4)) %>%
  gather(e_indicator, prop, c2_has1:c3_has4) %>%
  separate(e_indicator, c("class", "indicator"), "_")

test_seq_counts_cl <- t1 %>%
  mutate(rate_clust=gsub("c", "", values)) %>%
  summarise(seq=paste(rate_clust, collapse="")) %>% 
  group_by(ID, pop) %>% 
  do(data.frame(ID=.$ID, pop=.$pop, subseq=seqs, x=str_count(.$seq, seqs))) %>%
  group_by(pop, subseq) %>%
  summarise(n=sum(x)) %>%
  rowwise() %>%
  mutate(subseq=paste(unlist(strsplit(subseq, "")), collapse="_")) %>%
  separate(subseq, c("s1", "s2"), "_", remove=FALSE) %>%
  mutate(s1=as.numeric(s1), 
         s2=as.numeric(s2), 
         smin=pmin(s1,s2), 
         sord=as.numeric(paste0(smin, pmax(s1,s2))),
         has1=ifelse(s1==1|s2==1, TRUE, FALSE),
         has2=ifelse(s1==2|s2==2, TRUE, FALSE),
         has3=ifelse(s1==3|s2==3, TRUE, FALSE),
         has4=ifelse(s1==4|s2==4, TRUE, FALSE)) %>% #head
  gather(indicator, value, has1:has4) %>% #head
  dplyr::filter(value==TRUE & n>0) %>%
  group_by(pop, indicator, sord) %>%
  summarise(n=sum(n)) %>%
  group_by(pop, indicator) %>%
  mutate(tot=sum(n), prop=n/tot) %>% 
  mutate(class=ifelse(indicator=="has4", gsub("4", "", sord),
                      ifelse(indicator=="has3", gsub("3", "", sord),
                             ifelse(indicator=="has2", gsub("2", "", sord), gsub("1", "", sord))))) %>% 
  mutate(class=paste0("c", class)) %>% #head
  left_join(rc_counts, by=c("pop", "class", "indicator")) 

test_seq_counts_cl %>%
  rowwise() %>%
  do(tidy(prop.test(c(.$n), c(.$tot), p=.$prop.y)))


#-----------------------------------------------------------------------------
# plot histograms of cluster run frequencies
#-----------------------------------------------------------------------------
t1 %>%
  # dplyr::filter(rate_clust=="c1") %>%
  dplyr::filter(lengths<30) %>%
  ggplot(aes(x=lengths))+
  geom_histogram(binwidth=1)+
  facet_wrap(~values, scales="free")+
  # scale_x_log10()+
  theme_bw()

#-----------------------------------------------------------------------------
# test if relative rate per interval is associated with interval length 
#-----------------------------------------------------------------------------
test_sites2 <- test_sites %>% 
  dplyr::filter(D2>0 & D2<5e6) %>% 
  dplyr::filter(COUNT>0) %>%
  mutate(rate=COUNT/D2, rate2=rate^2, COUNT3=COUNT^3, D2m=pmin(D2, D2n)) %>% 
  # dplyr::filter(rate<0.05 & rate>0 & COUNT>0) %>% 
  # dplyr::filter(rate_clust!="c1") %>% 
  group_by(ID, pop) %>% 
  # mutate(z_score_rate = (rate - mean(rate)) / sd(rate)) %>% 
  do(augment(lm(log(.$D2)~.$rate+.$rate2)))
# do(augment(glm.nb(.$COUNT~.$D2m)))

test_sites2_summary <- test_sites %>% 
  dplyr::filter(D2>0 & D2<5e6) %>% 
  mutate(rate=COUNT/D2, rate2=rate^2, D2m=pmin(D2, D2n)) %>% 
  dplyr::filter(rate<0.05 & rate>0 & COUNT>0) %>% 
  # dplyr::filter(rate_clust!="c1") %>% 
  group_by(ID, pop) %>% 
  # mutate(z_score_rate = (rate - mean(rate)) / sd(rate)) %>% 
  do(tidy(lm(log(.$D2)~.$rate+.$rate2)))

msg <- "runtime"
pbPost("note", "test_sites2 lm finished", msg)

test_sites2 %>%
  # sample_n(5000) %>%
  dplyr::filter(ID=="NWD102555") %>%
  # dplyr::filter(.fitted>0 & .fitted<11) %>%
  dplyr::filter(..rate<0.025) %>%
  ggplot(aes(y=..rate, x=exp(log...D2.)))+
  geom_point(alpha=0.1)+
  geom_smooth(method="lm")+
  geom_vline(xintercept=2e4)+
  scale_x_log10()+
  scale_y_log10()+
  annotation_logticks()+
  xlab("distance")+
  ylab("relative mutation rate")

ggsave("ERV_mutation_hotspots/figs/NWD102555.dist.vs.rate.png", width=12, height=6)

test_sites2 %>% dplyr::filter(exp(.fitted)<20000) %>% nrow
nrow(test_sites2)

test_log <- function(x) {log(pexp(x, rate = 12000/3e9))}
test <- function(x) {pexp(x, rate = 12000/3e9)}
# test2 <- function(x) {pnorm(x, rate = 12000/3e9)}

# cdf
test_sites %>%
  dplyr::filter(pop=="EUR") %>%
  dplyr::filter(D2<1e6) %>%
  # dplyr::filter(rate_clust!="c1") %>%
  ggplot(.)+
  stat_ecdf(aes(x=D2), geom = "step", col="purple")+
  stat_ecdf(data=sim_sites, aes(x=dist), geom = "step", col="blue")+
  stat_ecdf(data=sim_sites_rc, aes(x=dist), geom = "step", col="grey30")+
  stat_ecdf(data=exp_df, aes(x=D2), geom = "step", col="green")+
  stat_function(fun=test, size=1, col="red")+
  scale_x_log10(limits=c(1,2e6))+
  xlab("Inter-singleton distance")+
  ylab("Cumulative probability")+
  annotation_logticks()+
  theme_bw()+
  theme()
# ggsave("ERV_mutation_hotspots/figs/compare_cdf.png", width=12, height=6)