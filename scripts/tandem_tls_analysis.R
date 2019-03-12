#-----------------------------------------------------------------------------
# test differences in tandem mutation spectra between pops 
#-----------------------------------------------------------------------------
test_sites %>% 
  ungroup() %>% 
  dplyr::filter(cl_ID!="UC") %>% 
  dplyr::filter(cl_LEN==cl_NUM) %>% 
  group_by(pop, TYPE) %>% 
  summarise(n=n()) %>% 
  mutate(tot=sum(n), prop=n/tot) %>% 
  group_by(TYPE) %>% 
  do(tidy(prop.test(.$n,.$tot))) %>% 
  data.frame

# mutation spectra of tandem mutations with 2,3,4,5+ sequential mutations
idsites %>% 
  dplyr::filter(cl_ID!="UC") %>% 
  dplyr::filter(cl_NUM==cl_LEN) %>% 
  ungroup() %>% 
  arrange(cl_ID) %>% 
  mutate(cl_NUM=ifelse(cl_NUM>=5, "5+", cl_NUM)) %>% 
  # group_by(pop, cl_NUM, TYPE) %>% 
  group_by(cl_NUM, TYPE) %>% 
  summarise(n=n()) %>% 
  # ungroup() %>%
  # mutate(cl_NUM2 = paste0(cl_NUM, " (", sum(n), ")")) %>% #head
  mutate(cl_NUM2=cl_NUM) %>%
  # group_by(pop, cl_NUM2) %>% 
  group_by(cl_NUM2) %>% 
  mutate(prop=n/sum(n)) %>% 
  ggplot(aes(x=TYPE, y=prop, fill=cl_NUM2))+
  geom_col(position="dodge")+
  scale_fill_manual("# bases mutated", values=brewer_pal(palette="YlOrRd")(6)[3:6])+
  # facet_wrap(~pop)+
  theme_bw()

#-----------------------------------------------------------------------------
# tls counts
#-----------------------------------------------------------------------------

tls_ids <- idsites %>%
  dplyr::filter(cl_LEN <= 100 & cl_NUM==2 & cl_ID!="UC") %>% 
  # mutate(TYPE=paste0(REF, "_", ALT, ".", substr(MOTIF, 3, 5))) %>%
  mutate(TYPE=paste0(REF, "_", ALT)) %>%
  dplyr::select(ID, CHR, POS, TYPE, cl_ID, cl_LEN) %>%
  separate(cl_ID, c("P1", "P2"), ":", remove=F) %>% #head # data.frame %>% #ungroup() %>%
  mutate(MNM_POS=ifelse(POS==P1, "P1", "P2")) %>% #head # mutate(POS_CHK=paste0(POS, ":")) %>%
  # dplyr::mutate(MNM_POS=ifelse(POS==gsub(":*", cl_ID), "P1", "P2")) %>% head
  #group_by(cl_ID, cl_LEN) %>% 
  dplyr::select(-c(POS, P1, P2)) %>%
  unite(cl_ID_new, ID, CHR, cl_ID, cl_LEN) %>% #head
  group_by(cl_ID_new) %>%
  spread(MNM_POS, TYPE) #%>% #head

tls_counts <- tls_ids %>%
  dplyr::filter(cl_ID_new != "NWD707063_chr20_10000058:9999998_-59") %>% # remove line with bug in cluster ID
  separate(cl_ID_new, c("ID", "CHR", "cl_ID", "cl_LEN"), sep="_") %>%
  mutate(MNM_TYPE=paste0(P1, ".", P2)) %>% #head(50)
  group_by(cl_LEN, P1, P2, MNM_TYPE) %>%
  count %>% #head
  # ungroup() %>%
  group_by(cl_LEN) %>%
  mutate(prop=n/sum(n))

# hm.palette <- colorRampPalette(rev(brewer.pal(9, 'PuGr')), space='Lab')

tls_counts %>%
  group_by(cl_LEN) %>%
  dplyr::filter(P1 %in% c("C_T", "T_C", "G_A", "A_G")) %>%
  dplyr::filter(P2 %in% c("C_T", "T_C", "G_A", "A_G")) %>%
  mutate(prop=n/sum(n)) %>%
  ungroup() %>%
  mutate(cl_LEN=as.numeric(cl_LEN)) %>% #head
  dplyr::filter(cl_LEN<=26) %>%
  ggplot(aes(x=P1, y=P2, fill=prop))+
  geom_tile()+
  facet_wrap(~factor(cl_LEN))+
  scale_fill_gradientn(colours=brewer.pal(7, 'Oranges'))

tls_counts %>%
  group_by(cl_LEN) %>%
  dplyr::filter(!P1 %in% c("C_T", "T_C", "G_A", "A_G")) %>%
  dplyr::filter(!P2 %in% c("C_T", "T_C", "G_A", "A_G")) %>%
  mutate(prop=n/sum(n)) %>%
  ungroup() %>%
  mutate(cl_LEN=as.numeric(cl_LEN)) %>% #head
  dplyr::filter(cl_LEN<=26) %>%
  ggplot(aes(x=P1, y=P2, fill=prop))+
  geom_tile()+
  facet_wrap(~factor(cl_LEN))+
  scale_fill_gradientn(colours=brewer.pal(7, 'Oranges'))

# tls_counts %>%
#   ggplot(aes(x=as.numeric(cl_LEN), y=prop, colour=MNM_TYPE))+
#   geom_point()+
#   scale_x_log10()+
# facet_wrap(~P1)
# gather(field, value, -c(cl_ID, cl_LEN)) %>% head(8)
# group_by(cl_ID, cl_LEN) %>% spread(field, value) %>% head
# spread()
# gather(TYPE,-c(TYPE, cl_LEN)) %>% head
# spread(cl_ID, TYPE) %>% head