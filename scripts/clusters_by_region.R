#-----------------------------------------------------------------------------  
# Model effects of genomic features
#-----------------------------------------------------------------------------  

if(!exists("cluster_counts_chr")){
  cluster_counts_chr <- test_sites %>%
    # dplyr::filter(rate_clust != "c1" | (cl_NUM<4 & rate_clust=="c1")) %>%
    mutate(BIN=ceiling(POS/1e6)) %>% # add bin column
    group_by(CHR, BIN, pop, rate_clust, ID) %>%
    summarise(n=n()) %>% # count singletons per bin, cluster, & ID
    group_by(CHR, BIN, pop, rate_clust) %>%
    dplyr::filter(n<mean(n)+2*sd(n)) %>% # in each bin, remove individuals with n>2sd from mean
    group_by(CHR, BIN, pop, rate_clust) %>%
    summarise(n=sum(n)) %>% # sum over all remaining individuals
    # mutate(min_dist=pmin(D2, D2n)) %>%
    # summarise(n=median(min_dist)) %>%
    group_by(pop, rate_clust) %>%
    mutate(prop=n/sum(n))
}

if(!exists("cluster_counts_chr_spectra")){
  cluster_counts_chr_spectra <- test_sites %>%
    # dplyr::filter(rate_clust != "c1" | (cl_NUM<4 & rate_clust=="c1")) %>%
    mutate(BIN=ceiling(POS/1e6)) %>%
    group_by(CHR, BIN, pop, rate_clust, ID, TYPE) %>%
    summarise(n=n()) %>% #head
    group_by(CHR, BIN, pop, rate_clust, ID) %>%
    mutate(tot=sum(n)) %>%
    group_by(CHR, BIN, pop, rate_clust) %>%
    dplyr::filter(n<mean(n)+2*sd(n)) %>%
    group_by(CHR, BIN, pop, rate_clust, TYPE) %>%
    summarise(n=sum(n)) %>% # sum by type over all individuals
    # mutate(min_dist=pmin(D2, D2n)) %>%
    # summarise(n=median(min_dist)) %>%
    group_by(CHR, BIN, pop, rate_clust) %>%
    mutate(prop=n/sum(n)) # calculate spectra
}

# cluster_counts_chr_spectra %>%
#   dplyr::filter(rate_clust=="c1" & TYPE %in% c("A_T", "C_A")) %>%
#   ggplot(aes(x=BIN, y=prop, colour=TYPE, group=TYPE))+
#   # geom_line()+
#   geom_point()+
#   facet_grid(CHR~pop)
# 
# cluster_counts_chr_spectra %>%
#   # dplyr::filter(CHR %in% c(2,8,9,16) & rate_clust=="c2") %>%
#   dplyr::filter(CHR %in% c(1:11) & rate_clust=="c2") %>%
#   ggplot(aes(x=BIN, y=prop, colour=TYPE, group=TYPE))+
#   # geom_line()+
#   geom_point()+
#   facet_grid(CHR~pop)

refdatadir, " -> refdatadir, "

h3k4me1 <- read_tsv(paste0(refdatadir, "/H3K4me1.1mb.bed"), col_names=F) %>%
  dplyr::select(CHR=X1, BIN=X2, h3k4me1=X7) %>%
  mutate(BIN=ceiling(BIN/1e6)+1)

h3k4me3 <- read_tsv(paste0(refdatadir, "/H3K4me3.1mb.bed"), col_names=F) %>%
  dplyr::select(CHR=X1, BIN=X2, h3k4me3=X7) %>%
  mutate(BIN=ceiling(BIN/1e6)+1)

h3k9ac <- read_tsv(paste0(refdatadir, "/H3K9ac.1mb.bed"), col_names=F) %>%
  dplyr::select(CHR=X1, BIN=X2, h3k9ac=X7) %>%
  mutate(BIN=ceiling(BIN/1e6)+1)

h3k9me3 <- read_tsv(paste0(refdatadir, "/H3K9me3.1mb.bed"), col_names=F) %>%
  dplyr::select(CHR=X1, BIN=X2, h3k9me3=X7) %>%
  mutate(BIN=ceiling(BIN/1e6)+1)

h3k27ac <- read_tsv(paste0(refdatadir, "/H3K27ac.1mb.bed"), col_names=F) %>%
  dplyr::select(CHR=X1, BIN=X2, h3k27ac=X7) %>%
  mutate(BIN=ceiling(BIN/1e6)+1)

h3k27me3 <- read_tsv(paste0(refdatadir, "/H3K27me3.1mb.bed"), col_names=F) %>%
  dplyr::select(CHR=X1, BIN=X2, h3k27me3=X7) %>%
  mutate(BIN=ceiling(BIN/1e6)+1)

h3k36me3 <- read_tsv(paste0(refdatadir, "/H3K36me3.1mb.bed"), col_names=F) %>%
  dplyr::select(CHR=X1, BIN=X2, h3k36me3=X7) %>%
  mutate(BIN=ceiling(BIN/1e6)+1)

rr <- read_tsv(paste0(refdatadir, "/rr.1mb.bed"), col_names=F) %>%
  dplyr::select(CHR=X1, BIN=X2, rr=X4) %>%
  dplyr::filter(rr != ".") %>%
  mutate(rr=as.numeric(rr)) %>%
  mutate(BIN=ceiling(BIN/1e6)+1)

cpgi <- read_tsv(paste0(refdatadir, "/cpgi.1mb.bed"), col_names=F) %>%
  dplyr::select(CHR=X1, BIN=X2, cpgi=X7) %>%
  mutate(BIN=ceiling(BIN/1e6)+1)

lamin <- read_tsv(paste0(refdatadir, "/lamin.1mb.bed"), col_names=F) %>%
  dplyr::select(CHR=X1, BIN=X2, lamin=X7) %>%
  mutate(BIN=ceiling(BIN/1e6)+1)

exon <- read_tsv(paste0(refdatadir, "/exons.1mb.bed"), col_names=F) %>%
  dplyr::select(CHR=X1, BIN=X2, exon=X7) %>%
  mutate(BIN=ceiling(BIN/1e6)+1)

dhs <- read_tsv(paste0(refdatadir, "/DHS.1mb.bed"), col_names=F) %>%
  dplyr::select(CHR=X1, BIN=X2, dhs=X7) %>%
  mutate(BIN=ceiling(BIN/1e6)+1)

features <- plyr::join_all(list(h3k4me1, 
                                h3k4me3, 
                                h3k9ac, 
                                h3k9me3, 
                                h3k27ac, 
                                h3k27me3, 
                                h3k36me3,
                                rr,
                                cpgi,
                                lamin,
                                exon,
                                dhs), by=c("CHR", "BIN"), type="full") %>%
  mutate(rr=replace_na(rr, 0))

cluster_features <- cluster_counts_chr %>%
  mutate(CHR=paste0(CHR)) %>%
  left_join(features, by=c("CHR", "BIN"))

cluster_features %>% 
  ungroup() %>%
  nest(-c(rate_clust, pop)) %>% 
  mutate(fit=purrr::map(data, ~ glm.nb(n~h3k4me1+h3k4me3+h3k9ac+h3k9me3+h3k27ac+h3k27me3+h3k36me3+rr+cpgi+lamin+exon+dhs, data=.x)), tidied=map(fit,tidy)) %>%
  unnest(tidied) %>%  #head
  data.frame %>% #dplyr::filter(p.value<0.05) %>%
  mutate(sig=ifelse(p.value<0.05, "*", " ")) %>%
  # mutate(estimate=exp(estimate)/10) %>%
  dplyr::filter(term != "(Intercept)") %>%
  ggplot(aes(x=term, y=estimate, fill=term))+
  geom_col()+
  geom_hline(yintercept=0)+
  geom_text(aes(label=sig, y=estimate*1.1))+
  xlab("Feature")+
  ylab("Beta")+
  # facet_grid(pop~rate_clust)+
  # facet_wrap(~rate_clust, ncol=1, strip.position="right")+
  facet_grid(rate_clust~pop)+
  theme_bw()+
  theme(axis.text.x=element_text(size=12, angle=45, hjust=1),
        legend.position="none",
        strip.text=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=16))#+
# ylab("Change in # singletons per 10% increase in feature")
ggsave(paste0(projdir, "/figs/feature_estimates.png"), width=12, height=8)

scor <- function(x){
  corrr::correlate(x, method="spearman")
}

#-----------------------------------------------------------------------------  
# c1 HMM
#-----------------------------------------------------------------------------  
mod <- depmix(response = prop ~ 1, 
              data = cluster_counts_chr[cluster_counts_chr$rate_clust=="c1",], 
              nstates = 3)

fm <- fit(mod)

cc1e <- cbind(data.frame(cluster_counts_chr[cluster_counts_chr$rate_clust=="c1",]), fm@posterior)
# cc1e %>% group_by(pop, state) %>% count()

cc1e %>%
  dplyr::filter(CHR %in% c("chr3", "chr4")) %>%
  mutate(CHR=paste0(CHR)) %>%
  dplyr::filter(prop<0.005) %>%
  mutate(state=factor(recode(state, "1"="cold", "2"="neutral", "3"="hot"), levels=c("hot", "neutral", "cold"))) %>%
  # dplyr::filter(pop == "AFR") %>%
  ggplot(aes(x=BIN, y=prop, fill=pop, alpha=factor(state)))+
  # geom_point(size=2, alpha=0.6)+
  geom_bar(stat="identity")+
  scale_alpha_manual("State", values=c(1,0.6,0.4))+
  # scale_fill_viridis("Ancestry", discrete=TRUE, begin=0, end=0.8)+
  scale_fill_manual("Ancestry", values=c("#FF7F00", "#33A02C"))+
  facet_grid(CHR~pop)+
  # facet_grid(CHR~rate_clust)+
  # scale_y_log10()+
  ylab("Proportion of \n component 1 singletons")+
  xlab("Physical position (Mb)")+
  guides(fill="none")+
  theme_classic()+
  # theme(legend.position="bottom")
  # ylab("Proportion of singletons")+
  # facet_grid(pop~rate_clust)+
  # facet_wrap(~rate_clust, ncol=1, strip.position="right")+
  theme(strip.text=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=16), 
        axis.text.x=element_text(size=12),
        # axis.ticks.x=element_blank(),
        # strip.text.y = element_text(angle = -90),
        legend.title=element_text(size=16),
        legend.text=element_text(size=12),
        legend.position="bottom")
ggsave(paste0(projdir, "/figs/chr3-4.component1.1mb.hmm.png"), width=12, height=4)

cc1e %>%
  # dplyr::filter(CHR %in% c(3:4)) %>%
  # mutate(CHR=factor(CHR, levels=1:22)) %>%
  # mutate(CHR=paste0("chr", CHR)) %>%
  mutate(CHR=factor(CHR, levels=paste0("chr", 1:22))) %>%
  dplyr::filter(prop<0.005) %>%
  mutate(state=factor(recode(state, "1"="cold", "2"="hot", "3"="neutral"), levels=c("hot", "neutral", "cold"))) %>%
  # dplyr::filter(pop == "AFR") %>%
  ggplot(aes(x=BIN, y=prop, fill=pop, alpha=factor(state)))+
  # geom_point(size=2, alpha=0.6)+
  geom_bar(stat="identity")+
  scale_alpha_manual("State", values=c(1,0.6,0.4))+
  # scale_fill_viridis("Ancestry", discrete=TRUE, begin=0, end=0.8)+
  scale_fill_manual("Ancestry", values=c("#FF7F00", "#E31A1C", "#33A02C"))+
  facet_grid(CHR~pop)+
  # facet_grid(CHR~rate_clust)+
  # scale_y_log10()+
  ylab("Proportion of \n component 1 singletons")+
  xlab("Physical position (Mb)")+
  guides(fill="none")+
  theme_classic()+
  # theme(legend.position="bottom")
  # ylab("Proportion of singletons")+
  # facet_grid(pop~rate_clust)+
  # facet_wrap(~rate_clust, ncol=1, strip.position="right")+
  theme(#strip.text=element_text(size=16),
    axis.title.y=element_text(size=16),
    axis.text.y=element_text(size=10),
    axis.title.x=element_text(size=16), 
    axis.text.x=element_text(size=12),
    strip.text.x=element_text(size=16),
    strip.text.y = element_text(size=16, angle = 0),
    # axis.ticks.x=element_blank(),
    legend.title=element_text(size=16),
    legend.text=element_text(size=12),
    legend.position="bottom")
# ggsave("ERV_mutation_hotspots/figs/component1.1mb.hmm.png", width=12, height=16)
ggsave(paste0(projdir, "/figs/component1.1mb.hmm.png"), width=12, height=12)

#-----------------------------------------------------------------------------  
# c2 HMM
#-----------------------------------------------------------------------------  
mod <- depmix(response = prop ~ pop, 
              data = cluster_counts_chr[cluster_counts_chr$rate_clust=="c2",], 
              # family = exponential(),
              nstates = 3)

fm <- fit(mod)

cc2e <- cbind(data.frame(cluster_counts_chr[cluster_counts_chr$rate_clust=="c2",]), fm@posterior)
cc2e %>% group_by(pop, state) %>% count()



c2_known <-  paste0("chr", c(2,8,9,16))

cc2e %>%
  group_by(pop) %>%
  mutate(rk=ntile(prop, 100)) %>%
  mutate(state=ifelse(rk>95, 1, 3)) %>%
  dplyr::filter(CHR %in% c2_known) %>%
  # mutate(CHR=paste0("chr", CHR)) %>%
  mutate(CHR=factor(CHR, levels=c2_known)) %>%
  mutate(state=factor(recode(state, "1"="hot", "2"="neutral", "3"="cold"), levels=c("hot", "neutral", "cold"))) %>%
  # dplyr::filter(pop == "AFR") %>%
  ggplot(aes(x=BIN, y=prop, fill=pop, alpha=factor(state)))+
  # geom_point(size=2, alpha=0.6)+
  geom_bar(stat="identity")+
  # scale_alpha_manual(values=c(1,0.4,0.2))+
  scale_alpha_manual("State", values=c(1,0.6,0.4))+
  # scale_fill_viridis("Ancestry", discrete=TRUE, begin=0, end=0.8)+
  scale_fill_manual("Ancestry", values=c("#FF7F00", "#E31A1C", "#33A02C"))+
  facet_grid(CHR~pop)+
  # facet_grid(CHR~rate_clust)+
  # scale_y_log10()+
  ylab("Proportion of \n component 2 singletons")+
  xlab("Physical position (Mb)")+
  guides(fill="none")+
  theme_classic()+
  # theme(legend.position="bottom")
  # ylab("Proportion of singletons")+
  # facet_grid(pop~rate_clust)+
  # facet_wrap(~rate_clust, ncol=1, strip.position="right")+
  theme(strip.text=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=16), 
        axis.text.x=element_text(size=12),
        strip.text.y = element_text(size=16, angle = 0),
        # axis.ticks.x=element_blank(),
        legend.title=element_text(size=16),
        legend.text=element_text(size=12),
        legend.position="bottom")
# ggsave("ERV_mutation_hotspots/figs/chr3-4.component1.1mb.hmm.png", width=12, height=4)
# ggsave("ERV_mutation_hotspots/figs/chr8-16.component2.2mb.hmm.png", width=12, height=12)
ggsave(paste0(projdir, "/figs/component2.chrs.1mb.hmm.png"), width=12, height=4)

cc2e %>%
  group_by(pop) %>%
  mutate(rk=ntile(prop, 100)) %>%
  mutate(state=ifelse(rk>95, 1, 3)) %>%
  # dplyr::filter(CHR %in% c(2,8,9,16)) %>%
  # mutate(CHR=paste0("chr", CHR)) %>%
  mutate(CHR=factor(CHR, levels=paste0("chr", 1:22))) %>%
  # mutate(CHR=factor(CHR, levels=c("chr2", "chr8", "chr9", "chr16"))) %>%
  mutate(state=factor(recode(state, "1"="hot", "2"="neutral", "3"="cold"), levels=c("hot", "neutral", "cold"))) %>%
  # dplyr::filter(pop == "AFR") %>%
  ggplot(aes(x=BIN, y=prop, fill=pop, alpha=factor(state)))+
  # geom_point(size=2, alpha=0.6)+
  geom_bar(stat="identity")+
  # scale_alpha_manual(values=c(1,0.4,0.2))+
  scale_alpha_manual("State", values=c(1,0.6,0.4))+
  # scale_fill_viridis("Ancestry", discrete=TRUE, begin=0, end=0.8)+
  scale_fill_manual("Ancestry", values=c("#FF7F00", "#E31A1C", "#33A02C"))+
  facet_grid(CHR~pop)+
  # facet_grid(CHR~rate_clust)+
  # scale_y_log10()+
  ylab("Proportion of \n component 2 singletons")+
  xlab("Physical position (Mb)")+
  guides(fill="none")+
  theme_classic()+
  # theme(legend.position="bottom")
  # ylab("Proportion of singletons")+
  # facet_grid(pop~rate_clust)+
  # facet_wrap(~rate_clust, ncol=1, strip.position="right")+
  theme(strip.text=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=16), 
        axis.text.x=element_text(size=12),
        strip.text.y = element_text(size=16, angle = 0),
        # axis.ticks.x=element_blank(),
        legend.title=element_text(size=16),
        legend.text=element_text(size=12),
        legend.position="bottom")
# ggsave("ERV_mutation_hotspots/figs/chr3-4.component1.1mb.hmm.png", width=12, height=4)
# ggsave("ERV_mutation_hotspots/figs/chr8-16.component2.2mb.hmm.png", width=12, height=12)
ggsave(paste0(projdir, "/figs/component2.1mb.hmm.png"), width=12, height=8)

#-----------------------------------------------------------------------------
# c3 HMM
#-----------------------------------------------------------------------------
mod <- depmix(response = prop ~ pop, 
              data = cluster_counts_chr[cluster_counts_chr$rate_clust=="c3",], 
              nstates = 3)

fm <- fit(mod)

cc3e <- cbind(data.frame(cluster_counts_chr[cluster_counts_chr$rate_clust=="c3",]), fm@posterior)
# cc3e %>% group_by(pop, state) %>% count()

cc3e %>%
  # dplyr::filter(CHR %in% c(2,8,9,16)) %>%
  # mutate(CHR=paste0("chr", CHR)) %>%
  mutate(CHR=factor(CHR, levels=paste0("chr", 1:22))) %>%
  # mutate(CHR=factor(CHR, levels=c("chr2", "chr8", "chr9", "chr16"))) %>%
  mutate(state=factor(recode(state, "1"="cold", "3"="neutral", "2"="hot"), levels=c("hot", "neutral", "cold"))) %>%
  # dplyr::filter(pop == "AFR") %>%
  ggplot(aes(x=BIN, y=prop, fill=pop, alpha=factor(state)))+
  # geom_point(size=2, alpha=0.6)+
  geom_bar(stat="identity")+
  # scale_alpha_manual(values=c(1,0.4,0.2))+
  scale_alpha_manual("State", values=c(1,0.6,0.4))+
  # scale_fill_viridis("Ancestry", discrete=TRUE, begin=0, end=0.8)+
  scale_fill_manual("Ancestry", values=c("#FF7F00", "#E31A1C", "#33A02C"))+
  facet_grid(CHR~pop)+
  # facet_grid(CHR~rate_clust)+
  # scale_y_log10()+
  ylab("Proportion of \n component 3 singletons")+
  xlab("Physical position (Mb)")+
  guides(fill="none")+
  theme_classic()+
  # theme(legend.position="bottom")
  # ylab("Proportion of singletons")+
  # facet_grid(pop~rate_clust)+
  # facet_wrap(~rate_clust, ncol=1, strip.position="right")+
  theme(strip.text=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=16), 
        axis.text.x=element_text(size=12),
        strip.text.y = element_text(size=16, angle = 0),
        # axis.ticks.x=element_blank(),
        legend.title=element_text(size=16),
        legend.text=element_text(size=12),
        legend.position="bottom")
# ggsave("ERV_mutation_hotspots/figs/chr3-4.component1.1mb.hmm.png", width=12, height=4)
# ggsave("ERV_mutation_hotspots/figs/chr8-16.component2.2mb.hmm.png", width=12, height=12)
ggsave(paste0(projdir, "/figs/component3.1mb.hmm.png"), width=12, height=12)

#-----------------------------------------------------------------------------
# c4 HMM
#-----------------------------------------------------------------------------
mod <- depmix(response = prop ~ pop, 
              data = cluster_counts_chr[cluster_counts_chr$rate_clust=="c4",], 
              nstates = 3)

fm <- fit(mod)

cc4e <- cbind(data.frame(cluster_counts_chr[cluster_counts_chr$rate_clust=="c4",]), fm@posterior)
cc4e %>% group_by(pop, state) %>% count()

cc4e %>%
  # dplyr::filter(CHR %in% c(2,8,9,16)) %>%
  # mutate(CHR=paste0("chr", CHR)) %>%
  mutate(CHR=factor(CHR, levels=paste0("chr", 1:22))) %>%
  # mutate(CHR=factor(CHR, levels=c("chr2", "chr8", "chr9", "chr16"))) %>%
  mutate(state=factor(recode(state, "3"="hot", "2"="neutral", "1"="cold"), levels=c("hot", "neutral", "cold"))) %>%
  # dplyr::filter(pop == "AFR") %>%
  ggplot(aes(x=BIN, y=prop, fill=pop, alpha=factor(state)))+
  # geom_point(size=2, alpha=0.6)+
  geom_bar(stat="identity")+
  # scale_alpha_manual(values=c(1,0.4,0.2))+
  scale_alpha_manual("State", values=c(1,0.6,0.4))+
  # scale_fill_viridis("Ancestry", discrete=TRUE, begin=0, end=0.8)+
  scale_fill_manual("Ancestry", values=c("#FF7F00", "#E31A1C", "#33A02C"))+
  facet_grid(CHR~pop)+
  # facet_grid(CHR~rate_clust)+
  # scale_y_log10()+
  ylab("Proportion of \n component 4 singletons")+
  xlab("Physical position (Mb)")+
  guides(fill="none")+
  theme_classic()+
  # theme(legend.position="bottom")
  # ylab("Proportion of singletons")+
  # facet_grid(pop~rate_clust)+
  # facet_wrap(~rate_clust, ncol=1, strip.position="right")+
  theme(strip.text=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=16), 
        axis.text.x=element_text(size=12),
        strip.text.y = element_text(size=16, angle = 0),
        # axis.ticks.x=element_blank(),
        legend.title=element_text(size=16),
        legend.text=element_text(size=12),
        legend.position="bottom")
# ggsave("ERV_mutation_hotspots/figs/chr3-4.component1.1mb.hmm.png", width=12, height=4)
# ggsave("ERV_mutation_hotspots/figs/chr8-16.component2.2mb.hmm.png", width=12, height=12)
ggsave(paste0(projdir, "/figs/component4.1mb.hmm.png"), width=12, height=12)

# test for differences
cc1e %>% 
  mutate(state=recode(state, "1"="3", "2"="2", "3"="1")) %>% 
  mutate(state=as.numeric(state)) %>%
  bind_rows(cc2e, cc3e, cc4e) %>% 
  mutate(state=factor(recode(state, "1"="hot", "2"="neutral", "3"="cold"), levels=c("hot", "neutral", "cold"))) %>%
  # mutate(CHR=paste0("chr", CHR)) %>%
  left_join(features, by=c("CHR", "BIN")) %>% 
  group_by(pop, rate_clust) %>%
  mutate(rk=ntile(prop, 100)) %>%
  dplyr::filter(rk>90 | rk<10) %>%
  mutate(state2=ifelse(rk>90, "hot", "cold")) %>%
  # dplyr::filter(prop>mean(prop)+2*sd(prop) | prop<mean(prop)-2*sd(prop)) %>% group_by(pop, rate_clust, state) %>% count
  dplyr::filter(state!="neutral") %>% 
  gather(feature, value, h3k4me1:dhs) %>% 
  group_by(pop, rate_clust, feature) %>% #head
  do(tidy(t.test(value~state2, data=.))) %>%
  mutate(sig=ifelse(p.value<0.05/11, TRUE, FALSE)) %>% #head
  dplyr::filter(sig==TRUE) %>% dplyr::filter(rate_clust=="c2")
# gather(state, estimate, estimate1:estimate2) %>%

