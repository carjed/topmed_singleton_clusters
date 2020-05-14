#-----------------------------------------------------------------------------
# run after selecting fullids from full file
#-----------------------------------------------------------------------------
# anc_ur <- anc %>% dplyr::filter(ID %in% fullids)
#
# afr_1k <- anc_ur %>%
#   dplyr::filter(AFR>0.85) %>%
#   arrange(desc(AFR)) %>%
#   head(1000) %>%
#   dplyr::select(ID)
#
# write_tsv(afr_1k, "/mnt/norbert/data/topmed/ancestry/afr_1k.txt", col_names=FALSE)
#
# eur_1k <- anc_ur %>%
#   dplyr::filter(EUR>0.85) %>%
#   arrange(desc(EUR)) %>%
#   head(1000) %>%
#   dplyr::select(ID)
#
# write_tsv(eur_1k, "/mnt/norbert/data/topmed/ancestry/eur_1k.txt", col_names=FALSE)

#-----------------------------------------------------------------------------
# read per-chromsome TOPMed data into single df--eats a lot of memory
#-----------------------------------------------------------------------------
# data_path <- "/mnt/norbert/data/topmed/singletons"   # path to the data
# sing_files <- dir(data_path, pattern = "*.sort.txt") # get file names
# # sing_files <- dir(data_path, pattern = "chr[1-9]") # get file names

# sing_files <- data.frame(name = sing_files) %>%
#   mutate(chr=getChrNum(name)) %>%
#   arrange(chr)
#
# sites <- sing_files$name[1:4] %>%
#   # read in all the files, appending the path before the filename
#   map(~ read_tsv(file.path(data_path, .))) %>%
#   reduce(rbind)

#-----------------------------------------------------------------------------
# TOPMed data
# - subset to 1000 individuals for more tractable computation
# - add column for 1Mb window
#-----------------------------------------------------------------------------

# sites <- read_tsv("/mnt/norbert/data/topmed/singletons/freeze3.singletons.sort.txt")
#
# testids <- head(unique(sites$ID), 1000)
# fullids <- unique(sites$ID)
# # length(testids)
#
# idsites <- sites %>%
#   dplyr::filter(ID %in% testids) %>%
#   mutate(BIN=paste0(CHR, ".", ceiling(POS/1e6)))

test_sites <- idsites %>%
  left_join(popids, by="ID") %>%
  left_join(sample_rates, by="ID") %>% #head
  left_join(avg_rates, by="pop") %>% #head
  # head() %>%
  # dplyr::filter(ID=="NWD100314") %>%
  group_by(ID) %>%
  dplyr::filter(D2>0) %>%
  mutate(D2n=lead(D2, n = 1L, default=1e6), Dmin=pmin(D2,D2n)) %>% #head
  # mutate(D2n=lead(D2, n = 1L, default=1e6)) %>% #head
  ungroup() %>% #head
  # group_by(CHR, POS) %>%
  mutate(rate_clust=pmap_chr(list(Dmin, p1, p2, p3, p4), assignCluster2)) %>% #,
  # rate_clust_m=map(assignCluster2(D2, D2n, c(p1a, p2a, p3a, p4a)))) %>% #head(20)
  # rowwise() %>%
  # mutate(rate_clust=assignCluster(D2, D2n, c(p1, p2, p3, p4)),
  #        rate_clust_m=assignCluster(D2, D2n, c(p1a, p2a, p3a, p4a))) %>% #head(20)
  # unnest(rate_clust) %>%
  ungroup() %>%
  dplyr::select(-c(p1:p4, p1a:p4a))
# mutate(rate_clust=assignCluster(D2, D2n, ID)) %>% head(20)

# exp_fits2 <- merge(exp_fits, id_counts, by="ID")
# exp_fits2_anc <- merge(exp_fits2, anc, by="ID") %>%
#   # dplyr::filter(n==4) %>%
#   gather(pop, pct_anc, AFR:ME) %>% arrange(ID) %>% #head
#   group_by(pop) %>%
#   mutate(ntile=ntile(pct_anc, 4)) %>%
#   dplyr::filter(pct_anc>0.5) %>%
#   dplyr::filter(pop %in% c("AFR", "EUR"))

# exp_fits2_anc %>%
# # rbind(exp_fits2_anc, exp_fits_br2) %>%
#   dplyr::filter(lambda>0.0001) %>%
#   # dplyr::filter(tot>14000) %>%
#   # dplyr::filter(tot>2000 & tot<10000) %>%
#   # ggplot(aes(x=rate, y=lambda, colour=param, shape=pop))+
#   ggplot(aes(x=rate, y=lambda, colour=pop, shape=param))+
#     scale_fill_gradient(low="grey50", high="white")+
#     # stat_density2d(geom = "raster", aes(fill = ..density..), contour = FALSE)+
#     # geom_point(aes(colour=.cluster), size=3, alpha=0.2)+
#     geom_point(size=3, alpha=0.6)+
#     # facet_wrap(~pop)+
#     # facet_grid(pop~ntile)+
#     # scale_colour_viridis(discrete=TRUE)+
#     scale_shape_manual(values=c(3, 1, 7, 6))+
#     # scale_shape(solid = FALSE)+
#     scale_x_log10(expand=c(0,0),
#                   breaks=c(0,1,10,100,1000,10000,100000,1e6),
#                   labels=c("0", "1", "10", "100", "1,000", "10,000", "100,000", "1M"))+
#     scale_y_log10(expand=c(0,0),
#                   limits=c(0.001,1.2),
#                   breaks=c(0.01, 0.1, 1))+
#     annotation_logticks()+
#     theme_bw()
#
# ggsave("ERV_mutation_hotspots/figs/exp.mix.rate.vs.prop.afr.eur.1k.png", width=12, height=6)

# exp_fits2_anc %>%
#   dplyr::filter(param %in% c("p2", "p3")) %>%
#   group_by(param, pop) %>%
#   summarise(cor=cor(rate, tot))

# Gaussian pdfs
# test_sites %>%
#   dplyr::filter(pop=="EUR") %>%
#   dplyr::filter(D2<3e6) %>%
#   # dplyr::filter(D2>100 & D2n>100) %>%
#   mutate(D2=log(D2)) %>%
#   ggplot(.)+
#   stat_function(fun=dnorm,
#                 args = list(mean = mean(log(exp_null)), sd = sd(log(exp_null))),
#                 size=1,
#                 linetype="dotted",
#                 # col="red",
#                 aes(colour="Expected (assuming uniform rate)"))+
#   geom_line(aes(x=D2, colour="Observed"),
#             # fill=NA,
#             # colour="blue",
#             size=1, stat="density")+
#   scale_x_continuous(limits=c(-1,15))+
#   scale_colour_manual("Model", values=c("red", "purple"))+
#   # scale_y_log10()+
#   xlab("log(inter-singleton distance)")+
#   ylab("density")+
#   # annotation_logticks()+
#   theme_bw()+
#   theme()
# ggsave("ERV_mutation_hotspots/figs/distance_distributions_gaussian.png", width=6, height=6)

# stat_ecdf(geom = "step", col="blue")+
# stat_function(fun=dnorm,
#               args = list(mean = mean(log(exp_p1)), sd = sd(log(exp_p1))),
#               size=1,
#               col="red")+
# stat_function(fun=dnorm,
#               args = list(mean = mean(log(exp_p2)), sd = sd(log(exp_p2))),
#               size=1,
#               col="red")+
# stat_function(fun=dnorm,
#               args = list(mean = mean(log(exp_p3)), sd = sd(log(exp_p3))),
#               size=1,
#               col="red")+
# stat_function(fun=dnorm,
#               args = list(mean = mean(log(exp_p4)), sd = sd(log(exp_p4))),
#               size=1,
#               col="red")+

# exp_fits2 %>%
#   as_tibble() %>%
#   dplyr::filter(lambda>0.0001) %>%
#   ggplot(aes(x=rate, y=lambda))+
#   scale_fill_gradient(low="grey50", high="white")+
#   stat_density2d(geom = "raster", aes(fill = ..density..), contour = FALSE)+
#   # geom_point(aes(colour=.cluster), size=3, alpha=0.2)+
#   geom_point(aes(colour=param), size=3, alpha=0.2)+
#   facet_wrap(~n)+
#   scale_colour_viridis(discrete=TRUE)+
#   scale_x_log10(expand=c(0,0),
#                 breaks=c(0,1,10,100,1000,10000,100000,1e6),
#                 labels=c("0", "1", "10", "100", "1,000", "10,000", "100,000", "1M"))+
#   scale_y_log10(expand=c(0,0),
#                 limits=c(0.001,1.2),
#                 breaks=c(0.01, 0.1, 1))+
#   annotation_logticks()+
#   theme_bw()
# ggsave("ERV_mutation_hotspots/figs/exp.mix.rate.vs.prop.windows.png", width=12, height=6)

#-----------------------------------------------------------------------------
# boxplots showing number of components vs number of singletons
#-----------------------------------------------------------------------------
# ggplot(exp_fits2, aes(x=factor(n), y=tot))+
#   geom_boxplot()+
#   xlab("number of components")+
#   ylab("singletons per sample (only chrs1-4)")+
#   theme_bw()
# ggsave("ERV_mutation_hotspots/figs/component.vs.nsing.boxplots.png", width=12, height=6)

#-----------------------------------------------------------------------------
# run k-means clustering
#-----------------------------------------------------------------------------
# kclusts3 <- exp_fits2 %>%
#   group_by(n) %>%
#   dplyr::select(lambda, rate) %>%
#   mutate(lambda=log10(lambda), rate=log10(rate)) %>%
#   do(kclust=kmeans(cbind(.$lambda,.$rate), mean(.$n)))
#
# exp_fits3c <- kclusts3 %>%
#   group_by(n) %>%
#   do(augment(.$kclust[[1]], exp_fits2[exp_fits2$n == .$n,]))
#
# exp_fits3_anc <- merge(exp_fits3c, anc, by="ID") %>%
#   dplyr::filter(n==4) %>%
#   gather(pop, pct_anc, AFR:ME) %>% arrange(ID) %>% #head
#   group_by(pop) %>%
#   mutate(ntile=ntile(pct_anc, 4)) %>%
#   dplyr::filter(pct_anc>0.5) %>%
#   dplyr::filter(pop %in% c("AFR", "EUR"))

# p <- exp_fits3_anc %>%
#   # rbind(exp_fits2_anc, exp_fits_br2) %>%
#   dplyr::filter(lambda>0.0001) %>%
#   # dplyr::filter(tot>2000 & tot<10000) %>%
#   # ggplot(aes(x=rate, y=lambda, colour=param, shape=pop))+
#   ggplot(aes(x=rate, y=lambda, colour=.cluster, shape=param))+
#     # scale_fill_gradient(low="grey50", high="white")+
#     # geom_point(aes(colour=.cluster), size=3, alpha=0.2)+
#     geom_point(size=3, alpha=0.6)+
#     stat_density2d(geom = "raster", aes(fill = ..density..), contour = FALSE)+
#     # facet_wrap(~pop)+
#     # facet_grid(pop~ntile)+
#     # scale_colour_viridis(discrete=TRUE)+
#     scale_shape_manual(values=c(3, 1, 7, 6))+
#     # scale_shape(solid = FALSE)+
#     scale_x_log10(expand=c(0,0),
#                   breaks=c(0,1,10,100,1000,10000,100000,1e6),
#                   labels=c("0", "1", "10", "100", "1,000", "10,000", "100,000", "1M"))+
#     scale_y_log10(expand=c(0,0),
#                   limits=c(0.001,1.2),
#                   breaks=c(0.01, 0.1, 1))+
#     annotation_logticks()+
#     theme_bw()

# p <- exp_fits2_anc %>%
#   # left_join(id_counts, by="ID") %>% head
#   dplyr::filter(!(pop == "EUR" & tot>15000)) %>%
#   # rbind(exp_fits2_anc, exp_fits_br2) %>%
#   dplyr::filter(lambda>0.0001) %>%
#   # dplyr::filter(tot>2000 & tot<10000) %>%
#   # ggplot(aes(x=rate, y=lambda, colour=param, shape=pop))+
#   # ggplot(aes(x=log10(rate), y=log10(lambda), colour=.cluster, shape=param))+
#   ggplot(aes(x=log10(rate), y=log10(lambda), group=pop, colour=pop, shape=param))+
#   scale_fill_gradient(low="grey20", high="white")+
#   stat_density2d(geom = "raster", aes(fill = ..density..), contour = FALSE)+
#   # geom_point(aes(colour=.cluster), size=3, alpha=0.2)+
#   # geom_point(size=3, alpha=0.6)+
#   # facet_wrap(~pop)+
#   # facet_grid(pop~ntile)+
#   scale_colour_viridis("Ancestry", discrete=TRUE, begin=0, end=0.8)+
#   scale_shape_manual("Mixture component", values=c(3, 1, 7, 6))+
#   scale_x_continuous(breaks=c(0:5), labels=10^c(0:5))+
#   scale_y_continuous(breaks=seq(-3, 0, 1), labels=10^seq(-3, 0, 1))+
#   xlab("1/rate")+
#   ylab("component contribution (lambda)")+
#   annotation_logticks()+
#   guides(colour=guide_legend(ncol=1),
#          shape=guide_legend(ncol=2))+
#   theme_bw()+
#   theme(axis.title.x=element_text(size=16),
#         axis.title.y=element_text(size=16),
#         axis.text.x=element_text(size=14),
#         axis.text.y=element_text(size=14),
#         # legend.position="bottom",
#         legend.position = c(0.9, 0.3),
#         legend.background = element_rect(fill="transparent", colour="black"))

# Gaussian pdfs comparing observed w/ 4-component mixtures
test_sites %>%
  # dplyr::filter(pop=="EUR") %>%
  group_by(pop) %>%
  sample_n(1e6) %>%
  dplyr::filter(D2<3e6) %>%
  # dplyr::filter(D2>100 & D2n>100) %>%
  mutate(D2=log(D2)) %>%
  ggplot(.)+
  geom_line(aes(x=D2, colour=pop, group=pop),
            # fill=NA,
            # colour="blue",
            size=1, stat="density")+
  geom_line(data=exp_df,
            aes(x=log(D2), colour="EUR"),
            # fill=NA,
            # colour="grey30",
            size=1, linetype="dashed", stat="density")+
  geom_line(data=exp_df_afr,
            aes(x=log(D2), colour="AFR"),
            # fill=NA,
            # colour="grey30",
            size=1, linetype="dashed", stat="density")+
  scale_x_continuous(limits=c(-1,15),
                     breaks=log(c(1,150,20000,300000)),
                     labels=c(1,150,20000,300000))+
  # scale_colour_manual("Model", values=c("darkgreen", "orange", "purple", "blue", "red"))+
  # scale_colour_manual("Model", values=c("darkgreen", "purple", "blue", "red"))+
  scale_colour_viridis("Ancestry", discrete=TRUE, begin=0, end=0.8)+
  # facet_wrap(~pop)+
  # scale_y_log10()+
  xlab("Inter-singleton distance")+
  ylab("density")+
  # annotation_logticks()+
  theme_bw()+
  theme(legend.position="bottom")
ggsave("ERV_mutation_hotspots/figs/distance_distributions_gaussian_mix.png", width=6, height=4)

# Gaussian pdfs
test_sites %>%
  # dplyr::filter(pop=="EUR") %>%
  group_by(pop) %>%
  sample_n(1e6) %>%
  dplyr::filter(D2<3e6) %>%
  # dplyr::filter(D2>100 & D2n>100) %>%
  mutate(D2=log(D2)) %>%
  ggplot(.)+
  geom_line(aes(x=D2, colour=pop, group=pop),
            # fill=NA,
            # colour="blue",
            size=1, stat="density")+
  stat_function(fun=dnorm,
                args = list(mean = mean(log(exp_null)), sd = sd(log(exp_null))),
                size=1,
                linetype="dotted",
                # col="red",
                aes(colour="EUR"))+  # expected (assuming exp(4.5e-06))
  stat_function(fun=dnorm,
                args = list(mean = mean(log(exp_null_afr)), sd = sd(log(exp_null_afr))),
                size=1,
                linetype="dotted",
                # col="red",
                aes(colour="AFR"))+ #  expected (assuming exp(6.3e-06))
  scale_x_continuous(limits=c(-1,15),
                     breaks=log(c(1,150,20000,300000)),
                     labels=c(1,150,20000,300000))+
  # scale_colour_manual("Model", values=c("darkgreen", "orange", "purple", "blue", "red"))+
  # scale_colour_manual("Model", values=c("darkgreen", "purple", "blue", "red"))+
  scale_colour_viridis("Ancestry", discrete=TRUE, begin=0, end=0.8)+
  # facet_wrap(~pop)+
  # scale_y_log10()+
  xlab("Inter-singleton distance")+
  ylab("density")+
  # annotation_logticks()+
  theme_bw()+
  theme(legend.position="bottom")
ggsave("ERV_mutation_hotspots/figs/distance_distributions_gaussian.png", width=6, height=6)


# Gaussian pdfs
test_sites %>%
  dplyr::filter(pop=="EUR") %>%
  dplyr::filter(D2<3e6) %>%
  # dplyr::filter(D2>100 & D2n>100) %>%
  mutate(D2=log(D2)) %>%
  ggplot(.)+
  stat_function(fun=dnorm,
                args = list(mean = mean(log(exp_null)), sd = sd(log(exp_null))),
                size=1,
                linetype="dotted",
                # col="red",
                aes(colour="1-component mixture"))+
  geom_line(data=exp_df,
            aes(x=log(D2), colour="4-component mixture"),
            # fill=NA,
            # colour="grey30",
            size=1, linetype="dotted", stat="density")+
  geom_line(data=sim_sites,
            aes(x=log(dist), colour="simulated (no RC)"),
            # fill=NA,
            # colour="grey30",
            size=1, linetype="dashed", stat="density")+
  geom_line(data=sim_sites_rc,
            aes(x=log(dist), colour="simulated (with RC)"),
            # fill=NA,
            # colour="grey30",
            size=1, linetype="dashed", stat="density")+
  # stat_density(aes(x=D2, colour="observed"),
  #              fill=NA,
  #              # colour="blue",
  #              size=1)+
  geom_line(aes(x=D2, colour="observed"),
            # fill=NA,
            # colour="blue",
            size=1, stat="density")+
  scale_x_continuous(limits=c(-1,15))+
  scale_colour_manual("Model", values=c("red", "green", "purple", "blue", "grey30"))+
  # scale_y_log10()+
  xlab("log(inter-singleton distance)")+
  ylab("density")+
  # annotation_logticks()+
  theme_bw()+
  theme()
ggsave("ERV_mutation_hotspots/figs/distance_distributions_gaussian.png", width=6, height=6)

#-----------------------------------------------------------------------------
# misc exploratory stuff
#-----------------------------------------------------------------------------
# exp_fits %>%
#   group_by(ID) %>%
#   dplyr::select(-lambda) %>%
#   spread(param, rate) %>%
#   dplyr::filter(abs(p1-p2)>1) %>%
#   gather(param, rate, p1:p3) %>%
#   arrange(ID, param)

# exp_fits_nparam <- exp_fits %>%
#   group_by(ID) %>%
#   summarise(n=n()) %>%
#   arrange(n)
#
# assignCluster2 <- function(dist, test_rates){
#   match <- which.max(lapply(test_rates, function(x) pexp(dist+1, 1/x)-pexp(dist-1,1/x)))
#   return(paste0("c", match))
# }

#-----------------------------------------------------------------------------
# misc plots
#-----------------------------------------------------------------------------
test_sites %>%
  dplyr::filter(D2<1e4 & D2>1e3) %>%
  # mutate(BIN=ceiling(POS/1e6)) %>%
  # group_by(CHR, BIN, pop) %>%
  # summarise(D2=mean(D2)) %>%
  ggplot(aes(x=D2, fill=pop))+
  # ggplot(aes(x=BIN, y=D2, colour=pop))+
  # geom_point()+
  # geom_bar(stat="identity", position="dodge")+
  # geom_density()+
  # scale_x_log10()+
  # scale_y_log10()+
  geom_histogram(aes(y = ..density..), position = 'identity', alpha=0.5, bins=50)+
  # facet_wrap(~CHR)+
  theme_bw()

test_sites %>% #data.frame %>% head
  mutate(MOTIF=substr(MOTIF,3,5)) %>%
  # dplyr::filter(!grepl("CG", motif)) %>%
  dplyr::filter(grepl("^A", TYPE)) %>%
  group_by(rate_clust, pop, TYPE, MOTIF) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  group_by(rate_clust, pop) %>%
  mutate(prop=n/sum(n)) %>% #head
  dplyr::filter(rate_clust=="c4") %>%
  # unite(subtype, TYPE:motif, sep="_") %>% #head
  ungroup() %>%
  dplyr::select(pop, TYPE, MOTIF, prop) %>% #head
  group_by(TYPE, MOTIF) %>%
  spread(pop, prop) %>% #head
  mutate(diff=EUR-AFR) %>%
  ungroup() %>%
  ggplot(aes(x=EUR, y=AFR))+
  geom_point(aes(colour=TYPE), alpha=0.5, size=4)+
  geom_smooth(method = "lm")+
  geom_abline(intercept=0, linetype="dashed")+
  scale_x_log10()+
  scale_y_log10()+
  xlab("fraction in European subsample")+
  ylab("fraction in African subsample")+
  theme_bw()

#-----------------------------------------------------------------------------
# count transitions/transversion per cluster/pop/1Mb window
#-----------------------------------------------------------------------------
tstv_counts_chr <- test_sites %>%
  # ungroup() %>%
  # slice(1:1e5) %>%
  mutate(BIN=ceiling(POS/1e6)) %>%
  group_by(CHR, BIN, pop, rate_clust, TYPE) %>%
  summarise(n=n()) %>%
  mutate(class=ifelse(TYPE %in% c("C_T", "A_G"), "TS", "TV")) %>%
  group_by(CHR, BIN, pop, rate_clust, class) %>%
  summarise(n=sum(n)) %>%
  group_by(CHR, BIN, pop, rate_clust) %>%
  mutate(prop=n/sum(n))

tstv_counts_chr %>%
  # dplyr::filter(rate_clust %in% c("c1", "c2")) %>%
  summarise(n=sum(n)) %>%
  # dplyr::filter(class=="TV") %>%
  dplyr::filter(n>30) %>%
  group_by(rate_clust, pop) %>%
  mutate(nnorm = n/mean(n)) %>%
  dplyr::filter(nnorm < 5) %>%
  # dplyr::filter(rate_clust %in% c("c1")) %>%
  # dplyr::filter(n<1000) %>%
  ggplot(aes(x=BIN, y=nnorm, colour=pop))+
  geom_point(size=2, alpha=0.6)+
  facet_grid(CHR~rate_clust, scales="free_y")+
  # facet_grid(CHR~rate_clust)+
  # scale_y_log10()+
  theme_bw()



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


# cluster_counts_chr %>%
#   dplyr::filter(rate_clust %in% c("c1")) %>%
#   dplyr::filter(prop<0.01) %>%
#   ggplot(aes(x=BIN, y=prop*100, fill=pop))+
#   # geom_point(size=2, alpha=0.6)+
#   geom_bar(stat="identity")+
#   scale_fill_viridis("Ancestry", discrete=TRUE, begin=0, end=0.8)+
#   facet_grid(CHR~pop)+
#   # facet_grid(CHR~rate_clust)+
#   # scale_y_log10()+
#   ylab("% of component 1 singletons")+
#   theme_classic()+
#   theme(legend.position="none")
# ggsave("ERV_mutation_hotspots/figs/singleton.component1.freqs.1mb.png", width=10, height=8)
#
# cluster_counts_chr %>%
#   dplyr::filter(rate_clust %in% c("c2")) %>%
#   ggplot(aes(x=BIN, y=prop*100, fill=pop))+
#     # geom_point(size=2, alpha=0.6)+
#     geom_bar(stat="identity")+
#     scale_fill_viridis("Ancestry", discrete=TRUE, begin=0, end=0.8)+
#     facet_grid(CHR~pop)+
#     # facet_grid(CHR~rate_clust)+
#     # scale_y_log10()+
#     ylab("% of component 2 singletons")+
#     theme_classic()+
#     theme(legend.position="none")
# ggsave("ERV_mutation_hotspots/figs/singleton.component2.freqs.1mb.png", width=10, height=8)

#
# NWD362844_fits <- idsites %>%
#   dplyr::filter(ID=="NWD362844") %>%
#   dplyr::filter(D2>0) %>%
#   # group_by(ID) %>%
#   do(fitExpMix(.$D2, scale=10, mincomp=3, maxcomp=4, iterate=FALSE)) %>%
#   as_tibble()
#
# test_sites %>%
#   group_by(pop, rate_clust) %>%
#   dplyr::summarise(n=n()) %>%
#   group_by(pop) %>% mutate(tot=sum(n), prop=n/tot)


#-----------------------------------------------------------------------------
# DEPRECATED
#-----------------------------------------------------------------------------
# # check mean values by state
# cc1e %>%
#   mutate(state=recode(state, "1"="3", "2"="2", "3"="1")) %>%
#   mutate(state=as.numeric(state)) %>%
#   bind_rows(cc2e, cc3e, cc4e) %>%
#   mutate(CHR=paste0("chr", CHR)) %>%
#   left_join(features, by=c("CHR", "BIN")) %>%
#   dplyr::filter(state!=2) %>%
#   gather(feature, value, h3k4me1:dhs) %>%
#   group_by(pop, rate_clust, feature, state) %>% #head
#   summarise(val=mean(value))
# # do(tidy(t.test(value~state, data=.))) %>%
# # mutate(sig=ifelse(p.value<0.05/11, TRUE, FALSE)) %>%
# # dplyr::filter(sig==TRUE)
#
# # boxplots
# cc1e %>%
#   mutate(state=recode(state, "1"="3", "2"="2", "3"="1")) %>%
#   mutate(state=as.numeric(state)) %>%
#   bind_rows(cc2e, cc3e, cc4e) %>%
#   # mutate(state=factor(recode(state, "1"="hot", "2"="neutral", "3"="cold"), levels=c("hot", "neutral", "cold"))) %>%
#   mutate(CHR=paste0("chr", CHR)) %>%
#   left_join(features, by=c("CHR", "BIN")) %>%
#   group_by(pop, rate_clust) %>%
#   # mutate(rk=ntile(prop, 100)) %>%
#   # dplyr::filter(rk>90 | rk<10) %>%
#   # mutate(state2=ifelse(rk>90, "hot", "cold")) %>% head
#   # dplyr::filter(prop>mean(prop)+2*sd(prop) | prop<mean(prop)-2*sd(prop)) %>% group_by(pop, rate_clust, state) %>% count
#   # dplyr::filter(state!="neutral") %>%
#   dplyr::filter(state!=2) %>%
#   gather(feature, value, h3k4me1:dhs) %>%
#   # mutate(state=recode(state, "1"="3", "2"="2", "3"="1")) %>%
#   # mutate(state=as.numeric(state)) %>%
#   # bind_rows(cc2e, cc3e, cc4e) %>%
#   # mutate(CHR=paste0("chr", CHR)) %>%
#   # left_join(features, by=c("CHR", "BIN")) %>%
#   # dplyr::filter(state!=2) %>%
#   # gather(feature, value, h3k4me1:dhs) %>%
#   group_by(pop, rate_clust, feature) %>% #head
#   ggplot(aes(x=feature, y=value, fill=factor(state)))+
#   geom_boxplot()+
#   facet_grid(rate_clust~pop)
# # do(tidy(t.test(value~state, data=.))) %>%
# # gather(state, estimate, estimate1:estimate2) %>%
#
# # logistic regression
# cc1e %>%
#   mutate(state=recode(state, "1"="3", "2"="2", "3"="1")) %>%
#   mutate(state=as.numeric(state)) %>%
#   bind_rows(cc2e, cc3e, cc4e) %>%
#   # mutate(state=factor(recode(state, "1"="hot", "2"="neutral", "3"="cold"), levels=c("hot", "neutral", "cold"))) %>%
#   mutate(CHR=paste0("chr", CHR)) %>%
#   left_join(features, by=c("CHR", "BIN")) %>%
#   # group_by(pop, rate_clust) %>%
#   # mutate(rk=ntile(prop, 100)) %>%
#   # dplyr::filter(rk>90 | rk<10) %>%
#   # mutate(state2=ifelse(rk>90, "hot", "cold")) %>% head
#   # dplyr::filter(prop>mean(prop)+2*sd(prop) | prop<mean(prop)-2*sd(prop)) %>% group_by(pop, rate_clust, state) %>% count
#   # dplyr::filter(state!="neutral") %>%
#   dplyr::filter(state!=2) %>%
#   mutate(state=recode(state, "1"="1", "3"="0")) %>%
#   mutate(state=as.numeric(state)) %>%
#   ungroup() %>%
#   nest(-c(rate_clust, pop)) %>%
#   mutate(fit=map(data, ~ glm(state~h3k4me1+h3k4me3+h3k9ac+h3k9me3+h3k27ac+h3k27me3+h3k36me3+cpgi+lamin+exon+dhs, data=.x, family="binomial")), tidied=map(fit,tidy)) %>%
#   # mutate(fit=map(data, ~ glm(prop~h3k4me1+h3k4me3+h3k9ac+h3k9me3+h3k27ac+h3k27me3+h3k36me3+cpgi+lamin+exon+dhs, data=.x)), tidied=map(fit,tidy)) %>%
#   unnest(tidied) %>% #head
#   data.frame %>% #dplyr::filter(p.value<0.05) %>%
#   mutate(sig=ifelse(p.value<0.05, "*", " ")) %>%
#   # mutate(estimate=exp(estimate)/10) %>%
#   dplyr::filter(term != "(Intercept)") %>%
#   ggplot(aes(x=term, y=estimate, fill=term))+
#   geom_col()+
#   geom_hline(yintercept=0)+
#   geom_text(aes(label=sig, y=estimate*1.1))+
#   xlab("Feature")+
#   ylab("Beta")+
#   # facet_grid(pop~rate_clust)+
#   # facet_wrap(~rate_clust, ncol=1, strip.position="right")+
#   facet_grid(rate_clust~pop)+
#   theme_bw()+
#   theme(axis.text.x=element_text(angle=45, hjust=1),
#         legend.position="none")#+

# cc2e %>%
#   dplyr::filter(CHR %in% c(8:16)) %>%
#   # dplyr::filter(pop == "AFR") %>%
#   ggplot(aes(x=BIN, y=n, fill=pop))+
#   # geom_point(size=2, alpha=0.6)+
#   geom_bar(stat="identity")+
#   # scale_alpha_manual(values=c(1,0.4,0.2))+
#   scale_fill_viridis("Ancestry", discrete=TRUE, begin=0, end=0.8)+
#   facet_grid(CHR~pop)+
#   # facet_grid(CHR~rate_clust)+
#   # scale_y_log10()+
#   ylab("# of component 2 singletons")+
#   xlab("Physical position (Mb)")+
#   theme_classic()+
#   theme(legend.position="none")
# ggsave("ERV_mutation_hotspots/figs/chr8-16.component2.2mb.png", width=12, height=12)
#
# cc2e %>% gather(class, prob, S1:S2) %>%
#   dplyr::filter(CHR %in% c(8:16)) %>%
#   # dplyr::filter(pop == "AFR") %>%
#   ggplot(aes(x=BIN, y=prob, colour=pop, alpha=class, group=class))+
#   # geom_point(size=2, alpha=0.6)+
#   geom_line(stat="identity")+
#   scale_alpha_manual(values=c(1,0.2))+
#   scale_colour_viridis("Ancestry", discrete=TRUE, begin=0, end=0.8)+
#   facet_grid(CHR~pop)+
#   # facet_grid(CHR~rate_clust)+
#   # scale_y_log10()+
#   # ylab("% of component 2 singletons")+
#   theme_classic()+
#   theme(legend.position="none")

# cluster_counts_chr %>%
#   dplyr::filter(rate_clust=="c2") %>%
#   group_by(pop) %>%
#   # dplyr::filter(prop<0.005) %>%
#   mutate(state=ifelse(n>mean(n)+2*sd(n), "hot", "cold")) %>%
#   # dplyr::filter(CHR %in% c(8:16)) %>%
#   # mutate(state=factor(recode(state, "1"="hot", "2"="neutral", "3"="cold"),
#   #                     levels=c("hot", "neutral", "cold"))) %>%
#   # mutate(state=factor(recode(state, "1"="cold", "2"="neutral", "3"="hot"),
#   #                     levels=c("cold", "neutral", "hot"))) %>%
#   # dplyr::filter(pop == "AFR") %>%
#   # ggplot(aes(x=BIN, y=n, fill=pop, alpha=factor(state)))+
#   ggplot(aes(x=BIN, y=prop*100, fill=pop, alpha=factor(state)))+
#     # geom_point(size=2, alpha=0.6)+
#     geom_bar(stat="identity")+
#     scale_alpha_manual(values=c(0.4, 1))+
#     scale_fill_viridis("Ancestry", discrete=TRUE, begin=0, end=0.8)+
#     facet_grid(CHR~pop)+
#     # facet_grid(CHR~rate_clust)+
#     # scale_y_log10()+
#     ylab("# of component 1 singletons")+
#     theme_classic()+
#     theme(legend.position="none")

# cc1e %>%
# cluster_counts_chr %>%
#   dplyr::filter(rate_clust=="c1") %>%
#   group_by(pop) %>%
#   # dplyr::filter(prop<0.005) %>%
#   mutate(state=ifelse(n>mean(n)+2*sd(n), "hot", "cold")) %>%
#   # dplyr::filter(CHR %in% c(8:16)) %>%
#   # mutate(state=factor(recode(state, "1"="hot", "2"="neutral", "3"="cold"), 
#   #                     levels=c("hot", "neutral", "cold"))) %>%
#   # mutate(state=factor(recode(state, "1"="cold", "2"="neutral", "3"="hot"), 
#   #                     levels=c("cold", "neutral", "hot"))) %>%
#   # dplyr::filter(pop == "AFR") %>%
#   ggplot(aes(x=BIN, y=n, fill=pop, alpha=factor(state)))+
#   # geom_point(size=2, alpha=0.6)+
#   geom_bar(stat="identity")+
#   scale_alpha_manual(values=c(0.4, 1))+
#   scale_fill_viridis("Ancestry", discrete=TRUE, begin=0, end=0.8)+
#   facet_grid(CHR~pop)+
#   # facet_grid(CHR~rate_clust)+
#   # scale_y_log10()+
#   ylab("# of component 1 singletons")+
#   theme_classic()+
#   theme(legend.position="none")
# ggsave("ERV_mutation_hotspots/figs/chr8-16.component1.2mb.hmm.png", width=12, height=12)

# cluster_features %>% 
#   # group_by(pop, CHR, BIN) %>%
#   # dplyr::summarise(n=sum(n)) %>%
#   dplyr::select(-CHR, -prop, -BIN) %>% 
#   # mutate(rate_clust="c1") %>%
#   group_by(pop, rate_clust) %>% 
#   nest() %>%
#   mutate(data = map(data, purrr::compose(scor))) %>% 
#   unnest() %>%  
#   dplyr::select(pop, rate_clust, feature=rowname, corr=n) %>%
#   dplyr::filter(feature != "n") %>%
#   ggplot(aes(x=feature, y=corr, fill=feature))+
#   geom_col()+
#   facet_grid(pop~rate_clust)+
#   theme_bw()+
#   theme(axis.text.x=element_text(angle=45, hjust=1))


# deprecated build37 ancestry for freeze 3 data
# anc <- read_delim("/mnt/norbert/data/topmed/ancestry/global_list_19815.txt", delim=" ", col_names=F)



  # deprecated freeze 3 data
  # exp_fits_fh <- "/mnt/norbert/data/topmed/fitmix/exp_fits.txt"

  # output from 4-component mode
  # exp_fits_fh <- "/mnt/norbert/data/topmed/fitmix/exp_fits_freeze5_min.txt"