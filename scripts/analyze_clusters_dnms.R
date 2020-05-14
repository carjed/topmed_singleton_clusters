#-----------------------------------------------------------------------------
# read decode DNMs
#-----------------------------------------------------------------------------

chr22_dnms_fh <- paste0(decodedatadir, "/decode_chr/anno/chr22.decode.dnms.sort.txt")

if(!file.exists(chr22_dnms_fh)){
  source("format_decode_dnms.R")
  # command(s) to process decode DNMs
}

data_path <- paste0(decodedatadir, "/decode_chr/anno/"   # path to the data
decode_dnm_files <- dir(data_path, pattern = "*.txt")

decode_dnm_files <- data.frame(name = decode_dnm_files) %>%
  mutate(chr=getChrNum(name)) %>%
  arrange(chr)

decode_dnms <- decode_dnm_files$name %>%
  # read in all the files, appending the path before the filename
  map(~ read_tsv(file.path(data_path, .), col_names=TRUE)) %>%
  reduce(rbind)

exp_fits_decode_dnm <- decode_dnms %>%
  dplyr::filter(D2>0) %>%
  do(fitExpMix(x=.$D2, scale=10, mincomp=3, maxcomp=4)) %>%
  as_tibble()

avg_rates_decode <- exp_fits_decode_dnm %>% 
  mutate(pop="EUR") %>%
  group_by(param, pop) %>%
  summarise(rate=median(rate)) %>% 
  ungroup() %>%
  mutate(param=paste0(param, "a")) %>% 
  group_by(pop) %>% 
  spread(param, rate)

dnms_decode_class <- decode_dnms %>%
  mutate(pop="EUR") %>%
  left_join(avg_rates_decode, by="pop") %>% #head
  group_by(ID) %>%
  dplyr::filter(D2>0) %>%
  mutate(D2n=lead(D2, n = 1L, default=1e7)) %>% #head
  ungroup() %>% #head
  rowwise() %>%
  mutate(rate_clust_m=assignCluster(D2, D2n, c(p1a, p2a, p3a))) %>% #head(20)
  ungroup() %>%
  dplyr::select(-c(p1a:p3a))

decode_dnms_v1 <- read_tsv(paste0(decodedatadir, "/decode_DNMs.tsv"), col_names=TRUE)
names(decode_dnms_v1)[1:5] <- c("CHR", "POS", "REF", "ALT", "ID")

decode_dnms_ages <- left_join(dnms_decode_class, decode_dnms_v1, by=c("CHR", "POS", "REF", "ALT", "ID")) %>%
  mutate(Dmin=pmin(D2, D2n)) %>%
  group_by(ID, rate_clust_m) %>%
  summarise(n=n(), F_age=median(Fathers_age_at_conception, na.rm=T), M_age=median(Mothers_age_at_conception, na.rm=T))



decode_dnms_ages_cl <- decode_dnms_ages %>% dplyr::filter(rate_clust_m=="c2")

mod <- glm.nb(n~F_age+M_age, data=decode_dnms_ages_cl)
summary(mod)
#-----------------------------------------------------------------------------
# Read freeze 5 DNMs
# New version filters to families, calls singletons in probands
#-----------------------------------------------------------------------------
data_path <- paste0(datadir, "/dnms/new_sorted/")   # path to the data
dnm_files <- dir(data_path, pattern = "*.txt") # get file names

dnm_files <- data.frame(name = dnm_files) %>%
  mutate(chr=getChrNum(name)) %>%
  arrange(chr)

dnms <- dnm_files$name %>%
  # read in all the files, appending the path before the filename
  map(~ read_tsv(file.path(data_path, .), col_names=TRUE)) %>%
  reduce(rbind)

dnms_anc <- merge(dnms, anc, by="ID")

# 57513 eur dnms
dnms_anc_eur <- dnms_anc %>% dplyr::filter(EUR>0.85)

mitchell_ids <- read_tsv(paste0(datadir, "/dnms/freeze5.mitchell.ids.txt"), col_names=F) %>%
  dplyr::select(ID=X1)

dnms_anc_eur_na <- dnms_anc_eur %>%
  dplyr::filter(!(ID %in% mitchell_ids$ID))

dnms_anc_eur_a <- dnms_anc_eur %>%
  dplyr::filter(ID %in% mitchell_ids$ID)

# 14256 afr dnms
dnms_anc_afr <- dnms_anc %>% dplyr::filter(AFR>0.85)
# 8435 eas dnms
dnms_anc_eas <- dnms_anc %>% dplyr::filter(EAS>0.85)

# eur/afr admixed?
# dnms_anc_adm <- dnms_anc %>% dplyr::filter(AFR>0.3 & EUR>0.3)

#-----------------------------------------------------------------------------
# Test for difference in mean mutation rate between EUR/AFR subsamples
#-----------------------------------------------------------------------------
bind_rows(dnms_anc_eur, dnms_anc_afr) %>%
  mutate(pop=ifelse(EUR>0.5, "EUR", "AFR")) %>%
  group_by(ID, pop) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  do(tidy(t.test(n~pop, data=.)))

#-----------------------------------------------------------------------------
# Plot distributions of #dnms in EUR/AFR
#-----------------------------------------------------------------------------
bind_rows(dnms_anc_eur, dnms_anc_afr) %>%
  mutate(pop=ifelse(EUR>0.5, "EUR", "AFR")) %>%
  group_by(ID, pop) %>%
  summarise(n=n()) %>%
  dplyr::filter(n<220) %>%
  ggplot(aes(x=n, group=pop, fill=pop))+
  geom_density(alpha=0.5)

dnm_spectra <- dnms %>% #data.frame %>% head
  mutate(motif=substr(MOTIF,3,5)) %>%
  group_by(TYPE, motif) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  mutate(prop=n/sum(n)) 

#-----------------------------------------------------------------------------
# Old (conservative) version--singletons only in probands
#-----------------------------------------------------------------------------
data_path_old <- paste0(datadir, "/dnms/old_sorted/"   # path to the data
dnm_files_old <- dir(data_path_old, pattern = "*.txt") # get file names

dnm_files_old <- data.frame(name = dnm_files_old) %>%
  mutate(chr=getChrNum(name)) %>%
  arrange(chr)

dnms_old <- dnm_files_old$name %>%
  # read in all the files, appending the path before the filename
  map(~ read_tsv(file.path(data_path_old, .), col_names=TRUE)) %>%
  reduce(rbind)

# compare mutation rates in eur/non-eur
dnms_old_anc <- merge(dnms_old, anc_5, by="ID")

# 57513 eur dnms
dnms_old_anc_eur <- dnms_old_anc %>% dplyr::filter(EUR>0.85)

# split for Amish vs non-Amish comparisons
dnms_old_anc_eur_na <- dnms_old_anc_eur %>%
  dplyr::filter(!(ID %in% mitchell_ids$ID))

dnms_old_anc_eur_a <- dnms_old_anc_eur %>%
  dplyr::filter(ID %in% mitchell_ids$ID)

# names(dnms_old) <- c("CHR", "POS", "REF", "ALT", "TYPE", "MOTIF", "ID")

dnm_spectra_old <- dnms_old %>% #data.frame %>% head
  mutate(motif=substr(MOTIF,3,5)) %>%
  group_by(TYPE, motif) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  mutate(prop=n/sum(n))

#-----------------------------------------------------------------------------
# Fit mixture distributions for DNMs
#-----------------------------------------------------------------------------
exp_fits_dnm <- dnms %>%
  dplyr::filter(D2>0) %>%
  do(fitExpMix(x=.$D2, scale=10, mincomp=3, maxcomp=4)) %>%
  as_tibble()

exp_fits_dnm_eur <- dnms_anc_eur %>%
  dplyr::filter(D2>0) %>%
  # mutate(D2n=lead(D2, n = 1L, default=1e7)) %>% #head
  mutate(D2n=lead(D2, n = 1L, default=1e7), Dmin=pmin(D2,D2n)) %>%
  do(fitExpMix(x=.$Dmin, scale=10, mincomp=3, maxcomp=3)) %>%
  as_tibble()

# non-amish
exp_fits_dnm_eur_na <- dnms_anc_eur_na %>%
  dplyr::filter(D2>0) %>%
  do(fitExpMix(x=.$D2, scale=10, mincomp=3, maxcomp=3)) %>%
  as_tibble()

# amish
exp_fits_dnm_eur_a <- dnms_anc_eur_a %>%
  dplyr::filter(D2>0) %>%
  do(fitExpMix(x=.$D2, scale=10, mincomp=3, maxcomp=3)) %>%
  as_tibble()

exp_fits_dnm_eur_sub <- dnms_anc_eur %>%
  dplyr::sample_n(nrow(dnms_anc_afr)) %>%
  dplyr::filter(D2>0) %>%
  do(fitExpMix(x=.$D2, scale=10, mincomp=3, maxcomp=3)) %>%
  as_tibble()

exp_fits_dnm_afr <- dnms_anc_afr %>%
  dplyr::filter(D2>0) %>%
  mutate(D2n=lead(D2, n = 1L, default=1e7), Dmin=pmin(D2,D2n)) %>%
  do(fitExpMix(x=.$Dmin, scale=10, mincomp=3, maxcomp=3)) %>%
  as_tibble()

# fails when we try per id
# exp_fits_dnm_afr_id <- dnms_anc_afr %>%
#   dplyr::filter(D2>0) %>%
#   group_by(ID) %>%
#   do(fitExpMix(x=.$D2, scale=10, mincomp=3, maxcomp=3)) %>%
#   as_tibble()

exp_fits_dnm_eas <- dnms_anc_eas %>%
  dplyr::filter(D2>0) %>%
  do(fitExpMix(x=.$D2, scale=10, mincomp=3, maxcomp=3)) %>%
  as_tibble()

exp_fits_dnm_old <-  dnms_old %>%
  dplyr::filter(D2>0) %>%
  do(fitExpMix(x=.$D2, scale=10, mincomp=3, maxcomp=4)) %>%
  as_tibble()

#-----------------------------------------------------------------------------
# assign each dnm to a cluster
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# AFR
#-----------------------------------------------------------------------------
avg_rates_afr <- exp_fits_dnm_afr %>% 
  mutate(pop="AFR") %>%
  group_by(param, pop) %>%
  summarise(rate=median(rate)) %>% 
  ungroup() %>%
  mutate(param=paste0(param, "a")) %>% 
  group_by(pop) %>% 
  spread(param, rate)

dnms_anc_afr_class <- dnms_anc_afr %>%
  mutate(pop="AFR") %>%
  left_join(avg_rates_afr, by="pop") %>% #head
  group_by(ID) %>%
  dplyr::filter(D2>0) %>%
  mutate(D2n=lead(D2, n = 1L, default=1e7)) %>% #head
  ungroup() %>% #head
  rowwise() %>%
  mutate(rate_clust_m=assignCluster(D2, D2n, c(p1a, p2a, p3a))) %>% #head(20)
  ungroup() %>%
  dplyr::select(-c(p1a:p3a))

#-----------------------------------------------------------------------------
# EAS
#-----------------------------------------------------------------------------
avg_rates_eas <- exp_fits_dnm_eas %>% 
  mutate(pop="EAS") %>%
  group_by(param, pop) %>%
  summarise(rate=median(rate)) %>% 
  ungroup() %>%
  mutate(param=paste0(param, "a")) %>% 
  group_by(pop) %>% 
  spread(param, rate)

dnms_anc_eas_class <- dnms_anc_eas %>%
  mutate(pop="EAS") %>%
  left_join(avg_rates_eas, by="pop") %>% #head
  group_by(ID) %>%
  dplyr::filter(D2>0) %>%
  mutate(D2n=lead(D2, n = 1L, default=1e7)) %>% #head
  ungroup() %>% #head
  rowwise() %>%
  mutate(rate_clust_m=assignCluster(D2, D2n, c(p1a, p2a, p3a))) %>% #head(20)
  ungroup() %>%
  dplyr::select(-c(p1a:p3a))

#-----------------------------------------------------------------------------
# EUR
#-----------------------------------------------------------------------------
avg_rates_eur <- exp_fits_dnm_eur %>% 
  mutate(pop="EUR") %>%
  group_by(param, pop) %>%
  summarise(rate=median(rate)) %>% 
  ungroup() %>%
  mutate(param=paste0(param, "a")) %>% 
  group_by(pop) %>% 
  spread(param, rate)

dnms_anc_eur_class <- dnms_anc_eur %>%
  mutate(pop="EUR") %>%
  left_join(avg_rates_eur, by="pop") %>% #head
  group_by(ID) %>%
  dplyr::filter(D2>0) %>%
  mutate(D2n=lead(D2, n = 1L, default=1e7)) %>% #head
  ungroup() %>% #head
  rowwise() %>%
  mutate(rate_clust_m=assignCluster(D2, D2n, c(p1a, p2a, p3a))) %>% #head(20)
  ungroup() %>%
  dplyr::select(-c(p1a:p3a))

dnms_anc_eur_class %>%
  mutate(cl20kb = ifelse(pmin(D2, D2n)<20000, "cDNM", "ncDNM")) %>%
  group_by(TYPE, cl20kb) %>%
  summarise(n=n()) %>%
  group_by(cl20kb) %>%
  mutate(prop=n/sum(n)) %>%
  # ggplot(aes(x=motif, y=prop, fill=rate_clust_m, alpha=pop))+
  ggplot(aes(x=TYPE, y=prop, fill=cl20kb))+
  geom_bar(stat="identity", position="dodge")+
  # scale_alpha_manual(values=c(0.3, 0.6,1))+
  facet_wrap(~TYPE, scales="free_x", nrow=1)+
  # facet_grid(TYPE~motif, scales="free_x")+
  # facet_grid(pop~rate_clust)+
  # facet_wrap(~rate_clust, ncol=1, strip.position="right")+
  theme_bw()+
  theme()

#-----------------------------------------------------------------------------
# Merge dnms with cluster assignment
#-----------------------------------------------------------------------------
# dnms_class <- bind_rows(dnms_anc_eur_class, dnms_anc_afr_class, dnms_anc_eas_class)
dnms_class <- bind_rows(dnms_anc_eur_class, dnms_anc_afr_class)

dnm_spectra <- dnms_class %>%
  mutate(motif=substr(MOTIF,3,5)) %>%
  group_by(rate_clust_m, pop, TYPE, motif) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  group_by(rate_clust_m, pop) %>%
  mutate(prop=n/sum(n)) %>%
  group_by(rate_clust_m, pop, TYPE) %>%
  summarise(n=sum(n)) %>%
  ungroup() %>%
  group_by(rate_clust_m, pop) %>%
  mutate(prop=n/sum(n)) 

# test for difference in spectra
dnm_spectra %>% 
  dplyr::select(-prop) %>% 
  # dplyr::filter(TYPE %in% c("A_C", "A_T", "C_G", "C_A")) %>%
  # dplyr::filter(TYPE %in% c("A_C", "A_T")) %>%
  group_by(rate_clust_m) %>% 
  spread(pop, n) %>% 
  do(tidy(chisq.test(cbind(.$AFR, .$EUR))))

#-----------------------------------------------------------------------------
# Plot dnm mutation spectra 
#-----------------------------------------------------------------------------
dnms_class %>% #dplyr::filter(TYPE=="C_G")
  mutate(motif=substr(MOTIF,3,5)) %>%
  # group_by(TYPE, motif, rate_clust_m, pop) %>%
  group_by(TYPE, rate_clust_m, pop) %>%
  summarise(n=n()) %>%
  group_by(rate_clust_m, pop) %>%
  mutate(prop=n/sum(n)) %>%
  mutate(TYPE=gsub("_", ">", TYPE)) %>%
  # ggplot(aes(x=motif, y=prop, fill=rate_clust_m, alpha=pop))+
  ggplot(aes(x=TYPE, y=prop, fill=rate_clust_m))+
  geom_bar(stat="identity", position="dodge")+
  scale_fill_manual(values=mpalette[c(1,2,3)])+
  # scale_alpha_manual(values=c(0.3, 1))+
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
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        legend.title=element_text(size=16),
        legend.text=element_text(size=12),
        legend.position="bottom")
ggsave(paste0(scriptdir, "/figs/dnm.spectra.new.png"), width=12, height=6)


## old version
# dnm_spectra %>% #dplyr::filter(TYPE=="C_G")
#   ggplot(aes(x=motif, y=prop, fill=TYPE))+
#   geom_bar(stat="identity", position="dodge")+
#   scale_alpha_manual(values=c(0.6,1))+
#   # facet_grid(rate_clust~TYPE, scales="free_x")+
#   # facet_grid(pop~rate_clust)+
#   facet_wrap(~TYPE, nrow=1, scales="free_x")+
#   theme_bw()+
#   theme(axis.text.x=element_text(angle=90))
# ggsave("ERV_mutation_hotspots/figs/dnm.spectra.new.png", width=12, height=6)
