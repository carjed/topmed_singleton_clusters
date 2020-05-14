rates <- read_tsv("/mnt/norbert/data/7bp_1000k_rates.txt")

# rates2 <- rates %>% 
#   mutate(Category2=gsub("cpg_", "", Category2)) %>% 
#   dplyr::select(Sequence, Category2, rel_prop) %>% 
#   spread(Category2, rel_prop, fill=0)

rates_fmt <- rates %>% 
  dplyr::select(Sequence, Category2, rel_prop) %>% 
  mutate(Category2=gsub("cpg_", "", Category2)) %>% 
  mutate(Category2=dplyr::recode(Category2, `AT_CG`="A_C", `AT_GC`="A_G", `AT_TA`="A_T", `GC_AT`="C_T", `GC_CG`="C_G", `GC_TA`="C_A")) %>% 
  mutate(subtype=paste0(Category2,  ".", substr(Sequence,1,7))) %>% 
  dplyr::select(subtype, rel_prop)

write_tsv(rates_fmt, "/mnt/norbert/data/rates.txt", col_names=F)
