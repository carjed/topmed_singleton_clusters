#-----------------------------------------------------------------------------
# This script is used to generate the lists of TOPMed samples to subset for
# multinucleotide mutation analysis
#-----------------------------------------------------------------------------

# read list of samples with duplicates removed
samples_55k <- read_tsv("/mnt/norbert/data/topmed/samples.phased.minDP0.remDuplicates.with_consent.no_controls.uniq_geno.txt", col_names = F)

# read data containing each sample's consent group
consent_groups <- read_tsv("/mnt/norbert/data/topmed/samples.notphased.minDP10.with_consent.no_controls.uniq_geno.txt", col_names = T)
names(consent_groups)[1] <- "X1"

# read global ancestry data
anc_55k <- read_delim("/mnt/norbert/data/topmed/freeze.5b.global.ancestry.txt", col_names = F, delim=" ")
anc_55k_header <- c("X1", "AFR", "CSA", "EAS", "EUR", "NAT", "OC", "ME")
names(anc_55k) <- anc_55k_header

# read list of related samples
relatives <- read_tsv("/mnt/norbert/data/topmed/relatives.txt", col_names = F)

# exclude related samples and samples from unconsented groups
samples_55k_consent <- left_join(samples_55k, consent_groups, "X1")
samples_55k_consent_anc <- left_join(samples_55k_consent, anc_55k, "X1")

samples_55k_consent_anc_unrel <- samples_55k_consent_anc %>%
  dplyr::filter(!(X1 %in% relatives$X1))

samples_55k_consent_anc_rel <- samples_55k_consent_anc %>%
  dplyr::filter(X1 %in% relatives$X1)

excluded_consent_groups <- c("DS-CS-RD", "DS-ASTHMA-IRB-MDS-RD", "DS-DHD-IRB-COL-NPU", "HMB-IRB-COL-NPU")

invalid_samples_rel <- samples_55k_consent_anc_rel %>%
  dplyr::filter(consent %in% excluded_consent_groups)

valid_samples <- samples_55k_consent_anc_unrel %>%
  dplyr::filter(!(consent %in% excluded_consent_groups)) #%>%
  # dplyr::filter(study %in% c("JHS", "GOLDN", "GeneSTAR", "HyperGEN", "MESA"))

valid_samples_rel <- samples_55k_consent_anc %>%
  dplyr::filter(!(consent %in% excluded_consent_groups))

probands <- read_tsv("/mnt/norbert/data/topmed/dnms/probands.txt")

valid_probands <- valid_samples_rel %>%
  dplyr::filter(X1 %in% probands$ID)

# read data containing PI info
sample_metrics_134k <- read_tsv("/mnt/norbert/data/topmed/TOPMed_metrics_2018-10-01-110352.tab.csv")
# names(consent_groups)[1] <- "X1"

sample_studies <- sample_metrics_134k %>% dplyr::select(X1=SAMPLE_ID, PI_NAME, STUDY2=STUDY)

permission_studies <- c("Arnett", "Barnes", "Burchard",
  "Correa", "Ellinor", "Goldn",
  "Mathias", "Mitchell", "Montgomery",
  "Ramachandran", "Redline", "Rotter", "Silverman")

# select IDs of 1000 European ancestry samples
samples_eur <- left_join(valid_samples, sample_studies, "X1") %>%
  dplyr::filter(PI_NAME %in% permission_studies) %>%
  # arrange(desc(EUR)) %>% head
  # dplyr::filter(STUDY2 %in% c("HyperGen", "Asthma_Afr", "Asthma", ))
  top_n(1000, EUR)

samples_eur %>% 
  dplyr::select(X1) %>% 
  write_tsv("/mnt/norbert/data/topmed/samples_55k_1k_eur.txt", col_names = F)

# select IDs of 1000 African ancestry samples
samples_afr <- left_join(valid_samples, sample_studies, "X1") %>%
  dplyr::filter(PI_NAME %in% permission_studies) %>%
  # arrange(desc(EUR)) %>% head
  # dplyr::filter(STUDY2 %in% c("HyperGen", "Asthma_Afr", "Asthma", ))
  top_n(1000, AFR)

samples_afr %>% 
  dplyr::select(X1) %>% 
  write_tsv("/mnt/norbert/data/topmed/samples_55k_1k_afr.txt", col_names = F)
