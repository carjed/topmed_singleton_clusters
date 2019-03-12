#-----------------------------------------------------------------------------
# plot mutation spectra grouped by pop & cluster
#-----------------------------------------------------------------------------
spectra_by_clust <- test_sites %>% #data.frame %>% head
  mutate(motif=substr(MOTIF,3,5)) %>%
  group_by(rate_clust, pop, TYPE, motif) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  group_by(rate_clust, pop) %>%
  mutate(prop=n/sum(n)) 

spectra_by_clust_1mer <- spectra_by_clust %>%
  group_by(rate_clust, pop, TYPE) %>%
  summarise(n=sum(n)) %>%
  ungroup() %>%
  group_by(rate_clust, pop) %>%
  mutate(prop=n/sum(n)) 

# test for difference in spectra
spectra_by_clust_1mer %>% 
  dplyr::select(-prop) %>% 
  # dplyr::filter(TYPE %in% c("A_C", "A_T", "C_G", "C_A")) %>%
  dplyr::filter(TYPE %in% c("A_C", "A_T")) %>%
  group_by(pop) %>% 
  spread(rate_clust, n) %>% 
  do(tidy(chisq.test(cbind(.$c3, .$c4))))

spectra_by_clust_1mer %>% 
  dplyr::select(-prop) %>% 
  # dplyr::filter(TYPE %in% c("A_C", "A_T", "C_G", "C_A")) %>%
  # dplyr::filter(TYPE %in% c("A_C", "A_T")) %>%
  group_by(rate_clust) %>% 
  spread(pop, n) %>% 
  do(tidy(chisq.test(cbind(.$AFR, .$EUR))))

devtools::install_github("awhstin/awtools")
library(awtools)

test_sites %>% #dplyr::filter(TYPE=="C_G")
  mutate(motif=substr(MOTIF,3,5)) %>%
  # group_by(TYPE, motif, rate_clust_m, pop) %>%
  group_by(TYPE, rate_clust, pop) %>%
  summarise(n=n()) %>%
  group_by(rate_clust, pop) %>%
  mutate(prop=n/sum(n)) %>%
  mutate(TYPE=gsub("_", ">", TYPE)) %>%
  ungroup() %>%
  mutate(rate_clust=gsub("c", "", rate_clust)) %>%
  # ggplot(aes(x=motif, y=prop, fill=rate_clust_m, alpha=pop))+
  ggplot(aes(x=TYPE, y=prop, fill=rate_clust))+
  geom_bar(stat="identity", position="dodge")+
  scale_fill_manual(values=mpalette[c(1,2,4,3)])+
  # scale_alpha_manual(values=c(0.3, 1))+
  # facet_wrap(~TYPE, scales="free_x", nrow=1)+
  facet_grid(pop~TYPE, scales="free_x")+
  guides(fill=guide_legend("Cluster class"))+
  ylab("Proportion of singletons")+
  # facet_grid(pop~rate_clust)+
  # facet_wrap(~rate_clust, ncol=1, strip.position="right")+
  theme_bw()+
  theme(strip.text=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.y=element_text(size=12),
        axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        legend.title=element_text(size=16),
        legend.text=element_text(size=12),
        legend.position="bottom")
ggsave("ERV_mutation_hotspots/figs/spectra.afr.eur.by.clust.png", width=12, height=4)

# 3-mer mutation spectra
spectra_by_clust %>% #dplyr::filter(TYPE=="C_G")
  ggplot(aes(x=motif, y=prop, fill=rate_clust, alpha=pop))+
  geom_bar(stat="identity", position="dodge")+
  scale_alpha_manual(values=c(0.6,1))+
  facet_grid(rate_clust~TYPE, scales="free_x")+
  # facet_grid(pop~rate_clust)+
  # facet_wrap(~rate_clust, ncol=1, strip.position="right")+
  theme_bw()+
  theme()
ggsave("ERV_mutation_hotspots/figs/spectra3.afr.eur.by.clust.png", width=12, height=6)

#-----------------------------------------------------------------------------
# plot change in spectra as inter-singleton distance increases
#-----------------------------------------------------------------------------
# spectra_by_dist <- left_join(idsites, popids, by="ID") %>% 
spectra_by_dist <- test_sites %>% 
  # group_by(ID) %>%
  # dplyr::filter(D2>0) %>%
  # mutate(D2n=lead(D2, n = 1L, default=1e6)) %>% #head
  # ungroup() %>% 
  # dplyr::filter(cl_ID!="UC") %>%
  # dplyr::select(-cl_LEN) %>%
  # dplyr::filter(rate_clust %in% c("c1", "c2")) %>% 
  # mutate(gp=ntile(cl_LEN, 50)) %>% 
  # mutate(cl_LEN=pmin(D2, D2n)) %>% #arrange(desc(cl_LEN))
  mutate(cl_LEN=Dmin) %>% #arrange(desc(cl_LEN))
  dplyr::filter(cl_LEN<100000) %>%
  mutate(cl_LEN=ifelse(cl_LEN>100 & cl_LEN<20000, ceiling(cl_LEN/100)*100, cl_LEN)) %>%
  # mutate(cl_LEN=ifelse(cl_LEN>1000 & cl_LEN<10000, ceiling(cl_LEN/100)*100, cl_LEN)) %>%
  mutate(cl_LEN=ifelse(cl_LEN>20000, ceiling(cl_LEN/1000)*1000, cl_LEN)) %>%
  group_by(cl_LEN, TYPE, pop) %>% 
  summarise(n=n()) %>%
  group_by(cl_LEN, pop) %>% 
  mutate(prop=n/sum(n), tot=sum(n), error = sqrt((prop * (1-prop))/n)) #%>%
# dplyr::filter(cl_LEN<=10000) %>%
# mutate(TYPE_d=paste0(TYPE, ceiling(cl_LEN/10))) %>%
# left_join(t1, by=c("TYPE", "cl_LEN"))

spectra_by_dist %>%
  mutate(TYPE=gsub("_", ">", TYPE)) %>%
  # dplyr::filter(cl_LEN<=100) %>%
  mutate(gp=ifelse(cl_LEN<=100, "1-100bp", ifelse(cl_LEN<20000, "101-20,000bp", ">20,000bp"))) %>%
  mutate(gp=factor(gp, levels=c("1-100bp", "101-20,000bp", ">20,000bp"))) %>%
  # mutate(breakpt=ifelse(is.na(breakpt), ""))
  ggplot(aes(x=cl_LEN, y=prop, color=TYPE, fill=TYPE, group=TYPE))+
  # scale_x_log10(breaks=sort(c(5*10^(0:3), 10^(0:3))))+
  # scale_x_log10(expand=c(0,0.02), breaks=sort(c(1:9, seq(10,100,by=10))))+
  # geom_point(aes(size=tot), alpha=0.5)+
  scale_x_log10(breaks=c(1,10,100,1000,10000,50000, 100000), labels=1/c(1,10,100,1000,10000,50000, 100000))+
  # scale_x_continuous(breaks=seq(0,500,25))+
  geom_line(alpha=0.5)+
  # geom_vline(xintercept=14, linetype="dashed")+
  geom_smooth(span=0.1, se=F)+
  geom_point(size=0.5, alpha=0.5)+
  # geom_point(aes(shape=breakpt), size=4)+
  # scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,0.4,by=0.05))+
  scale_shape_manual(values=c(16,0), guide=FALSE)+
  facet_grid(pop~gp, scales="free_x")+
  # facet_wrap(~gp, scales="free_x")+
  # geom_errorbar(aes(ymin = prop - 1.96*error, ymax = prop + 1.96*error, color=TYPE), width=0.05)+
  # geom_bar(stat="identity")+
  # geom_vline(data=t1, aes(xintercept=rn, colour=TYPE))+
  # geom_vline(data=t1, aes(xintercept=x))+
  # annotation_logticks(sides="b")+
  # xlab("distance to nearest singleton (bp)")+
  xlab("Intrinsic mutation rate")+
  ylab("Proportion of singletons")+
  # guides(colour=guide_legend("Mutation Type"))+
  theme_classic()+
  theme(panel.border=element_rect(colour="black", fill=NA),
        strip.text=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=16),
        axis.text.x=element_text(size=12),
        # axis.title.x=element_blank(), 
        # axis.text.x=element_blank(), 
        # axis.ticks.x=element_blank(),
        legend.title=element_text(size=16),
        legend.text=element_text(size=12),
        legend.position="bottom")
ggsave("ERV_mutation_hotspots/figs/spectra_by_dist.png", width=12, height=8)