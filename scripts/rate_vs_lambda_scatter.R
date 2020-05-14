#-----------------------------------------------------------------------------
# rate vs. lambda scatterplots
#-----------------------------------------------------------------------------

# plot how number of singletons per sample correlates with intrinsic mutation rates
exp_fits2_anc %>% 
  dplyr::filter(param %in% c("p1", "p2", "p3")) %>% 
  ggplot(aes(x=rate, y=tot, group=pop, colour=pop, shape=param))+
  geom_point(size=3, alpha=0.6)+
  geom_smooth(method="lm")+
  scale_colour_viridis("Ancestry", discrete=TRUE, begin=0, end=0.8)+
  # scale_x_log10()+
  scale_shape_manual("Mixture component", values=c(3, 1, 7))+
  # scale_x_continuous(breaks=c(0:5), labels=10^c(0:5))+
  facet_wrap(~param, scales="free")+
  # scale_y_continuous(breaks=seq(-3, 0, 1), labels=10^seq(-3, 0, 1))+
  xlab("1/rate")+
  ylab("Number of singletons")+
  # annotation_logticks()+
  guides(colour=guide_legend(ncol=1))+
  theme_bw()+
  theme(axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        # legend.position="bottom",
        legend.position = "bottom",
        legend.background = element_rect(fill="transparent", colour="black"))

p2_quartiles <- exp_fits2_anc %>%
  dplyr::filter(n==4) %>%
  # left_join(id_counts, by="ID") %>% head
  # dplyr::filter(!(pop == "EUR" & tot>15000)) %>%
  left_join(exp_fits2_anc %>% 
              dplyr::filter(param=="p1") %>% 
              group_by(pop) %>%
              # mutate(rank=ntile(tot, 100)) %>% 
              mutate(rank=ntile(tot, 4)) %>% 
              ungroup() %>%
              # dplyr::filter(rank>10) %>%
              dplyr::select(ID, rank), by="ID") %>% #head
  # group_by(pop) %>%
  # mutate(rank=ntile(tot, 100)) %>%
  # mutate(group=ifelse(rank>75, "top25", "bottom75")) %>%
  # mutate(group=ntile(tot, 4)) %>%
  # dplyr::filter(tot>12000) %>%
  # rbind(exp_fits2_anc, exp_fits_br2) %>%  
  dplyr::filter(lambda>0.0001) %>%
  mutate(param=gsub("p", "", param)) %>%
  # dplyr::filter(tot>2000 & tot<10000) %>%
  # ggplot(aes(x=rate, y=lambda, colour=param, shape=pop))+
  # ggplot(aes(x=log10(rate), y=log10(lambda), colour=.cluster, shape=param))+
  ggplot(aes(x=log10(rate), y=log10(lambda), group=pop, colour=pop, shape=param))+
  # scale_fill_gradient(low="grey50", high="white")+
  # stat_density2d(geom = "raster", aes(fill = ..density..), contour = FALSE)+
  # geom_point(aes(colour=.cluster), size=3, alpha=0.2)+
  geom_point(size=3, alpha=0.6)+
  # facet_wrap(~pop)+
  # facet_grid(pop~ntile)+
  # scale_colour_viridis("Ancestry", discrete=TRUE, begin=0, end=0.8)+
  scale_colour_manual("Ancestry", values=c("#FF7F00", "#E31A1C", "#33A02C"))+
  scale_shape_manual("Cluster class", values=c(3, 1, 7, 6, 4))+
  # scale_x_continuous(breaks=c(0:5), labels=10^c(0:-5))+
  scale_x_continuous(breaks=c(0:5), labels=10^c(0:5))+
  scale_y_continuous(breaks=seq(-3, 0, 1), labels=10^seq(-3, 0, 1))+
  facet_wrap(~rank, ncol=1)+
  # xlab("1/rate")+
  # xlab("Intrinsic mutation rate")+
  xlab("Inter-singleton distance (bp)")+
  ylab("component contribution (lambda)")+
  annotation_logticks()+
  guides(colour=guide_legend(ncol=1),
         shape=guide_legend(ncol=2))+
  theme_bw()+
  theme(axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        # legend.position="bottom",
        legend.position = c(0.15, 0.7),
        legend.background = element_rect(fill="transparent", colour="black"))

# save without marginal histograms
p2_quartiles
ggsave(paste0(projdir, "/figs/exp.mix.rate.vs.prop.afr.eur.1k.filter.png", width=12, height=12)

p2 <- exp_fits2_anc %>%
  dplyr::filter(n==4) %>%
  # left_join(id_counts, by="ID") %>% head
  # dplyr::filter(tot>15000) %>%
  # dplyr::filter(!(pop == "EUR" & tot>15000)) %>%
  left_join(exp_fits2_anc %>% 
              dplyr::filter(param=="p1") %>% 
              group_by(pop) %>%
              # mutate(rank=ntile(tot, 100)) %>% 
              mutate(rank=ntile(tot, 4)) %>% 
              ungroup() %>%
              # dplyr::filter(rank>10) %>%
              dplyr::select(ID, rank), by="ID") %>% #head
  # group_by(pop) %>%
  # mutate(rank=ntile(tot, 100)) %>%
  # mutate(group=ifelse(rank>75, "top25", "bottom75")) %>%
  # mutate(group=ntile(tot, 4)) %>%
  # dplyr::filter(tot>12000) %>%
  # rbind(exp_fits2_anc, exp_fits_br2) %>%  
  dplyr::filter(lambda>0.0001) %>%
  mutate(param=gsub("p", "", param)) %>%
  # dplyr::filter(tot>2000 & tot<10000) %>%
  # ggplot(aes(x=rate, y=lambda, colour=param, shape=pop))+
  # ggplot(aes(x=log10(rate), y=log10(lambda), colour=.cluster, shape=param))+
  ggplot(aes(x=log10(rate), y=log10(lambda), group=pop, colour=pop, shape=param))+
  # scale_fill_gradient(low="grey50", high="white")+
  # stat_density2d(geom = "raster", aes(fill = ..density..), contour = FALSE)+
  # geom_point(aes(colour=.cluster), size=3, alpha=0.2)+
  geom_point(size=2, alpha=0.6)+
  # facet_wrap(~pop)+
  # facet_grid(pop~ntile)+
  # scale_colour_viridis("Ancestry", discrete=TRUE, begin=0, end=0.8)+
  scale_colour_manual("Ancestry", values=c("#FF7F00", "#E31A1C", "#33A02C"))+
  scale_shape_manual("Cluster class", values=c(3, 1, 4, 0, 4))+
  # scale_x_continuous(breaks=c(0:5), labels=10^c(0:-5))+
  scale_x_continuous(breaks=c(0:5), labels=10^c(0:5))+
  scale_y_continuous(breaks=seq(-3, 0, 1), labels=10^seq(-3, 0, 1))+
  # facet_wrap(~rank, ncol=1)+
  # xlab("1/rate")+
  # xlab("Intrinsic mutation rate")+
  xlab("Inter-singleton distance (bp)")+
  ylab("component contribution (lambda)")+
  annotation_logticks()+
  guides(colour=guide_legend(ncol=1),
         shape=guide_legend(ncol=2))+
  theme_bw()+
  theme(axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        # legend.position="bottom",
        legend.position = c(0.15, 0.7),
        legend.background = element_rect(fill="transparent", colour="black"))

# save without marginal histograms
# p2
# ggsave("ERV_mutation_hotspots/figs/exp.mix.rate.vs.prop.afr.eur.1k.filter.png", width=12, height=12)

# save with marginal histograms
marg.hist.title <- paste0(projdir, "/figs/exp.mix.rate.vs.prop.afr.eur.1k.marg.hist.png")
png(marg.hist.title, width=12, height=6, units="in", res=300)
ggMarginal(p2, type="histogram", position="identity", groupFill=T, colour=NA,  alpha = 0.8, xparams = list(bins=100), yparams = list(bins=100))
dev.off()

#-----------------------------------------------------------------------------
# Generate GIF of scatterplot, adding first 50 individuals
#-----------------------------------------------------------------------------
p2a <- exp_fits2_anc %>%
  dplyr::filter(ID %in% unique(exp_fits2_anc$ID)[1:50]) %>%
  # left_join(id_counts, by="ID") %>% head
  dplyr::filter(!(pop == "EUR" & tot>15000)) %>%
  # rbind(exp_fits2_anc, exp_fits_br2) %>%  
  dplyr::filter(lambda>0.0001) %>%
  # dplyr::filter(tot>2000 & tot<10000) %>%
  # ggplot(aes(x=rate, y=lambda, colour=param, shape=pop))+
  # ggplot(aes(x=log10(rate), y=log10(lambda), colour=.cluster, shape=param))+
  ggplot(aes(x=log10(rate), y=log10(lambda), group=pop, colour=pop, shape=param))+
  # scale_fill_gradient(low="grey50", high="white")+
  # stat_density2d(geom = "raster", aes(fill = ..density..), contour = FALSE)+
  # geom_point(aes(colour=.cluster), size=3, alpha=0.2)+
  geom_point(size=3)+
  # facet_wrap(~pop)+
  # facet_grid(pop~ntile)+
  # scale_colour_viridis("Ancestry", discrete=TRUE, begin=0, end=0.8)+
  scale_colour_manual("Ancestry", values=c("#FF7F00", "#E31A1C", "#33A02C"))+
  scale_shape_manual("Mixture component", values=c(3, 1, 7, 6))+
  scale_x_continuous(breaks=c(0:5), labels=10^c(0:5))+
  scale_y_continuous(breaks=seq(-3, 0, 1), labels=10^seq(-3, 0, 1))+
  xlab("1/rate")+
  ylab("component contribution (lambda)")+
  annotation_logticks()+
  guides(colour=guide_legend(ncol=1),
         shape=guide_legend(ncol=2))+
  theme_bw()+
  theme(axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        # legend.position="bottom",
        legend.position = c(0.15, 0.7),
        legend.background = element_rect(fill="transparent", colour="black"))+
  transition_manual(ID, cumulative=TRUE)

anim_save(filename = paste0(projdir, "/figs/mixture_component_animation.gif"), animation = p2a)
