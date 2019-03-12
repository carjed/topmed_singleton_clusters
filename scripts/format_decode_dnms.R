# Read single-nucleotide DNMs 

decode_dnms <- read_tsv("/mnt/norbert/data/dnms/decode_DNMs.tsv") %>% 
  # dplyr::filter(is.na(Discordant_in_3_gen_or_mz_twins)) %>% #head
  dplyr::filter(str_length(Ref) == 1 & str_length(Alt) == 1) %>% #head(40)
  mutate(CAT=paste0(Ref, Alt)) %>%
  mutate(TYPE=CAT)

decode_dnms$TYPE[decode_dnms$CAT=="AC" | decode_dnms$CAT=="TG"] <- "A_C"
decode_dnms$TYPE[decode_dnms$CAT=="AG" | decode_dnms$CAT=="TC"] <- "A_G"
decode_dnms$TYPE[decode_dnms$CAT=="AT" | decode_dnms$CAT=="TA"] <- "A_T"
decode_dnms$TYPE[decode_dnms$CAT=="GA" | decode_dnms$CAT=="CT"] <- "C_T"
decode_dnms$TYPE[decode_dnms$CAT=="GC" | decode_dnms$CAT=="CG"] <- "C_G"
decode_dnms$TYPE[decode_dnms$CAT=="GT" | decode_dnms$CAT=="CA"] <- "C_A"

decode_dnms <- decode_dnms %>%
  mutate(MOTIF=substr(TYPE, 1, 1)) %>%
  dplyr::select(CHR=Chr, POS=Pos_hg38, REF=Ref, ALT=Alt, TYPE, MOTIF, ID=Proband_nr)

decode_dnms %>%
  group_by(CHR) %>%
  do(write_tsv(., paste0("/mnt/norbert/data/dnms/decode_chr/", unique(.$CHR), ".decode.dnms.txt"), col_names=F))

