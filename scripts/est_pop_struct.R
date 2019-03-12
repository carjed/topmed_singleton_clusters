#-----------------------------------------------------------------------------
# use singleton distance distribution to capture population structure
#-----------------------------------------------------------------------------

exp_fits_pca2 <- test_sites %>% 
  group_by(ID, rate_clust, TYPE) %>%
  count() %>%
  group_by(ID) %>%
  mutate(tot=sum(n), prop=n/tot) %>%
  mutate(param3=paste0(rate_clust, "_", TYPE)) %>%
  ungroup() %>%
  dplyr::select(ID, param3, prop) %>%
  spread(param3, prop) %>% #head
  left_join(exp_fits2_anc %>%
              ungroup() %>%
              dplyr::filter(n==4) %>%
              # dplyr::filter(param != "p1") %>%
              gather(param2, value, c(lambda, rate)) %>%
              mutate(param3=paste0(param, "_", param2)) %>%
              dplyr::select(ID, pct_anc, pop, param3, value) %>%
              spread(param3, value),
            by="ID") %>%
  # left_join(exp_fits2_anc %>% dplyr::filter(param=="p1") %>% dplyr::select(pop, pct_anc, ID), by="ID") %>%
  nest() %>%
  mutate(pca = map(data, ~prcomp(.x %>% dplyr::select(-c(pop, pct_anc, ID)), center = TRUE, scale = TRUE)), 
         pca_tidy = map2(pca, data, ~broom::augment(.x, data = .y))) #%>% unnest() %>% head

exp_fits_pca2[[3]][[1]] %>% head
  
# left_join(exp_fits2_anc %>% 
#             ungroup() %>%
#             dplyr::filter(n==4) %>%
#             # dplyr::filter(param != "p1") %>%
#             gather(param2, value, c(lambda, rate)) %>% 
#             mutate(param3=paste0(param, "_", param2)) %>% 
#             dplyr::select(ID, pct_anc, pop, param3, value) %>% 
#             spread(param3, value),
#           by="ID") %>%  

exp_fits_pca <- exp_fits2_anc %>% 
  ungroup() %>%
  dplyr::filter(n==4) %>%
  # dplyr::filter(param != "p1") %>%
  gather(param2, value, c(lambda, rate)) %>% 
  mutate(param3=paste0(param, "_", param2)) %>% 
  dplyr::select(ID, pct_anc, pop, param3, value) %>% 
  spread(param3, value) %>% #head
  nest() %>% 
  mutate(pca = map(data, ~prcomp(.x %>% dplyr::select(-c(pop, pct_anc, ID)), center = TRUE, scale = FALSE)), 
         pca_tidy = map2(pca, data, ~broom::augment(.x, data = .y))) #%>% unnest() %>% head


exp_fits_umap <- exp_fits_pca2[[3]][[1]] %>% 
  dplyr::select(.fittedPC1:.fittedPC8) %>%
  data.frame() %>%
  umap()

bind_cols(exp_fits_pca2[[3]][[1]], data.frame(exp_fits_umap$layout)) %>% 
  left_join(anc, by="ID") %>%
  ggplot(aes(x=X1, y=X2, colour=pop))+
  geom_point(alpha=0.5)+
  scale_colour_manual("Ancestry", values=c("#FF7F00", "#33A02C"))+
  theme_bw()

exp_fits_pca2[[3]][[1]] %>% 
  ggplot(aes(x=.fittedPC1, y=.fittedPC2, colour=pop, alpha=pct_anc))+
  geom_point()+
  theme_bw()
