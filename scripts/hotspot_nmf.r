setwd("./ERV_mutation_hotspots/")

source("http://bioconductor.org/biocLite.R")
biocLite("Biostrings")

library(Biostrings)
library(viridis)
library(heatmaply)
library(dendextend)

################################################################################
# Gene-wise analysis
################################################################################

spectra <- read.table("data/mut_sigs_genes/NMF_M_spectra.txt", header=T, stringsAsFactors=F)
spectra_rates <- read.table("data/mut_sigs_genes/NMF_M_spectra_rates.txt", header=T, stringsAsFactors=F)
sig_loads <- read.table("data/mut_sigs_genes/NMF_H_sig_loads.txt", header=T, stringsAsFactors=F)
sig_contribs <- read.table("data/mut_sigs_genes/NMF_W_sig_contribs.txt", header=T, stringsAsFactors=F)


sig_contribs2 <- sig_contribs[complete.cases(sig_contribs),] %>%
  # mutate(sumVar = rowSums(.[2:ncol(sig_contribs)])) %>%
  mutate(sumVar = rowSums(.[2:ncol(sig_contribs)])) %>%
  mutate_at(vars(starts_with("S")), funs(./sumVar)) %>%
  mutate(maxgp = apply(.[,2:ncol(sig_contribs)], 1, function(x)
    names(x)[which.max(x)])) %>%
  dplyr::select(-sumVar) %>%
  gather(signature, contribution, -c(ID,maxgp))


sig_loads_long <- sig_loads %>%
  mutate(sumVar = rowSums(.[2:ncol(sig_loads)])) %>%
  rowwise() %>%
  # mutate_each(funs(./sumVar), -Sig) %>%
  mutate_at(vars(contains(".")), funs(./sumVar)) %>%
  gather(subtype, loading, -c(Sig, sumVar)) %>%
  separate(subtype, c("category", "motif"), sep = "[.]")

ggplot(sig_loads_long, aes(x=motif, y=loading, fill=Sig))+
  geom_bar(stat="identity")+
  facet_grid(Sig~category, scales="free")+
  scale_fill_viridis(discrete=TRUE)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, vjust=0, size=10),
        axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24),
        strip.text=element_text(size=24),
        legend.position="none")
ggsave("figs/gene_wise_sigloads.png", width=12.5, height=7.3)

genes_list <- sig_contribs2 %>% 
  filter(maxgp=="S5") %>% 
  arrange(desc(contribution)) %>% 
  head(10) %>% 
  dplyr::select(ID) %>% 
  unlist

maxgp <- sig_contribs2 %>% 
  dplyr::filter(contribution > 0.4) %>%
  dplyr::select(ID, maxgp) 

# s5genes <- "PAX1, SOX1, FOXC1, GSX1, OTX1, MSX1, DMRTA2, HOXB3, HOXC4, HOXD3"
# s5genes_list <- strsplit(gsub(" ", "", s5genes), ',') %>% unlist

genes <- read.table("data/gencode.v19.names.bed", header=F, stringsAsFactors=F)
names(genes) <- c("CHR", "START", "END", "ID")
sc2 <- merge(sig_contribs, genes, by="ID")

# sc2 <- sc2[!anyDuplicated(sc2$ID),]
gene_ERV_counts <- spectra %>% gather(subtype, count, -ID) %>% group_by(ID) %>% summarise(n=sum(count))
gene_rates <- merge(gene_ERV_counts, genes, by="ID")
gene_rates <- gene_rates[!(duplicated(gene_rates$ID) | duplicated(gene_rates$ID, fromLast = TRUE)), ]

sc2 <- sc2[!(duplicated(sc2$ID) | duplicated(sc2$ID, fromLast = TRUE)), ]
sc2$LENGTH = sc2$END - sc2$START

# spectra3 <- spectra %>%
spectra3 <- merge(spectra, maxgp, by="ID") %>%
  # mutate(group=ifelse(ID %in% keep_ids$V1, "keep", "drop")) %>%
  group_by(maxgp) %>%
# spectra3 <- spectra %>%
  gather(subtype, count, -c(ID, maxgp)) %>%
  separate(subtype, c("category", "motif"), sep = "[.]") %>%
  # group_by(ID, category, motif) %>%
  group_by(maxgp, category, motif) %>%
  summarise(count=sum(count)) %>%
  group_by(maxgp) %>%
  mutate(prop=count/sum(count)) %>%
  ungroup() 


sig_contribs2 %>% 
  dplyr::filter(ID %in% c("MCPH1", "DEFA4", "DEFA5", "DEFB1", "CSMD1")) %>% 
  ggplot(aes(x=ID, y=contribution, fill=signature))+
    geom_bar(stat="identity", position="stack")+
  coord_flip()+
  scale_fill_viridis(discrete=TRUE)+
  # scale_x_discrete(expand=c(0,0))+
  # scale_y_discrete(labels=seq(0,1,0.1), breaks=seq(0,1,0.1))+
  xlab("Gene")+
  theme_bw()+
  theme(legend.position="bottom",
    axis.title.x=element_text(size=16),
    axis.title.y=element_text(size=16),
    strip.text.y=element_text(size=16),
    strip.text.x=element_text(size=16),
    axis.text.y=element_text(size=12),
    axis.text.x=element_text(size=12))

ggsave("figs/chr8p_genes_sig_contribs.png", width=10, height=7) 

top_sig_genes <- sig_contribs2 %>% 
  dplyr::filter(!grepl("RP|AC0|CTD|LINC", ID)) %>% 
  group_by(maxgp) %>% 
  arrange(desc(contribution)) %>% 
  slice(1L) %>% 
  ungroup() %>% 
  dplyr::select(ID) %>% unlist

spectra_rates %>%
  dplyr::filter(ID %in% top_sig_genes) %>%
  gather(Subtype, contribution, -ID) %>% 
  separate(Subtype, into=c("Type", "Motif"), sep="[.]") %>%
  ggplot(aes(x=Motif, y=contribution, fill=ID))+
    geom_bar(stat="identity")+
    facet_grid(ID~Type, scales="free")+
    scale_fill_brewer(palette = "Set1")+
    theme_bw()+
    theme(axis.text.x=element_text(angle=90, vjust=0, size=10),
          axis.text.y=element_text(size=20),
          axis.title.x=element_text(size=24),
          axis.title.y=element_text(size=24),
          strip.text=element_text(size=24),
          legend.position="none")
ggsave("figs/sample_spectra.png", width=14, height=7)

# spectra3 <- spectra3 %>%
#   filter(category %in% c("A_C", "A_T", "C_A", "C_G"))
# 
# s3 <- spectra3 %>% dplyr::select(-count) %>% spread(maxgp, prop)
# 
# heatmaply(cor(s3[,3:7]), margins = c(40, 40))

# ggplot(spectra3, aes(x=motif, y=prop, fill=ID, group=ID))+
ggplot(spectra3, aes(x=motif, y=prop, fill=maxgp, group=maxgp))+
  geom_bar(stat="identity", position="dodge")+
  facet_grid(maxgp~category, scales="free_x")+
  # facet_wrap(~category, ncol=6, scales="free_x")+
  ylab("Contribution")+
  scale_fill_brewer(palette = "Set1")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, vjust=0, size=20),
        axis.text.y=element_text(size=10),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        legend.title=element_blank(),
        legend.text=element_text(size=10),
        strip.text=element_text(size=12))

################################################################################
# 1Mb window analysis
################################################################################

sigloads <- read.table("data/sigloads.txt", header=T, stringsAsFactors=F)
bins_nmf_long <- read.table("data/bin_sig_prob.txt", header=T, stringsAsFactors=F)
ind_nmf_long <- read.table("data/ind_sig_prob.txt", header=T, stringsAsFactors=F)

bins_nmf_long2 <- bins_nmf_long %>%
  dplyr::select(CHR, BIN, window, X1:X5) %>%
  mutate(sumVar = rowSums(.[4:8])) %>%
  mutate_at(vars(starts_with("X")), funs(./sumVar)) %>%
  group_by(window) %>%
  slice(1L) %>%
  gather(cluster, prob_m, X1:X5)

mean_contrib <- bins_nmf_long2 %>% 
  group_by(cluster) %>% 
  summarise(mean_sig=mean(prob_m), sd_sig=sd(prob_m))

sig_peaks_bed <- merge(bins_nmf_long2, mean_contrib, by="cluster") %>% 
  rowwise() %>% 
  dplyr::filter(prob_m>(mean_sig+2*sd_sig)) %>% 
  mutate(START=BIN*1e6-1e6, END=BIN*1e6) %>% 
  dplyr::select(CHR, START, END, cluster)

# full_data$sites %>%
#   dplyr::filter(CHR==sig_peaks_bed[1,]$CHR & POS<sig_peaks_bed[1,]$END & POS>sig_peaks_bed[1,]$START) %>%
#   group_by(ID, Category) %>%
#   count() 

write.table(sig_peaks_bed, "data/sig_peaks.bed", col.names=F, row.names=F, sep="\t", quote=F)

################################################################################
# using NMF matrices generated by python script
################################################################################

nmf_M <- read.csv("/mnt/norbert/projects/ERV_mutation_hotspots/data/nmf_1Mb_M.csv")
names(nmf_M)[1] <- "ID"

nmf_M1 <- nmf_M %>% 
  gather(subtype, rate, -c(ID)) %>%
  separate(subtype, c("type", "motif"), sep="[.]")

orderedcats1a <- c("A_C", "A_G", "A_T",
                  "C_A", "C_G", "C_T")

orderedcats2a <- c("A>G", "A>C", "A>T",
                  "C>A", "C>G", "C>T")

myPaletteCat <- colorRampPalette(brewer.pal(12, "Paired"))

gp_cols <- myPaletteCat(12)[
  c(10,8,12,
    2,4,6,
    1,3,5)]


# rates7$Type <- ifelse(substr(rates7$Motif, 4, 5)=="CG",
#                       paste0("cpg_", rates7$Type), rates7$Type)

nmf_M1$Category2 <- plyr::mapvalues(nmf_M1$type, orderedcats1a, orderedcats2a)
nmf_M1$Category2 <- factor(nmf_M1$Category2, levels=orderedcats2a)

int_breaks <- function(x, n = 5) pretty(x, n)[pretty(x, n) %% 1 == 0]

ggplot(nmf_M1, aes(x=log(rate), y=motif, fill=type, colour=type, height=..density..))+
  geom_density_ridges(stat="binline", bins=50, draw_baseline = FALSE, alpha=0.8)+
  scale_x_continuous(breaks=int_breaks)+
  scale_fill_manual(values=gp_cols, drop=FALSE)+
  scale_colour_manual(values=gp_cols, drop=FALSE)+
  facet_wrap(~Category2, scales="free", nrow=1)+
  theme_bw()+
  theme(legend.position="none",
        # axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        strip.text.y=element_text(size=16),
        strip.text.x=element_text(size=16))
ggsave("/mnt/norbert/projects/ERV_mutation_hotspots/figs/rate_window_dist.png", width=12, height=5.25)

nmf_H <- read.csv("/mnt/norbert/projects/ERV_mutation_hotspots/data/nmf_1Mb_H.csv")
names(nmf_H)[1] <- "signature"

bg_rates <- read.csv("/mnt/norbert/projects/ERV_mutation_hotspots/data/bg_rates.csv")
names(bg_rates)[1] <- "subtype"

bg_rates <- bg_rates %>%
  separate(subtype, c("type", "motif"), sep="[.]")

sl1 <- nmf_H %>%
  mutate(sumVar = rowSums(.[2:97])) %>%
  mutate_at(vars(contains("_")), funs(./sumVar)) %>%
  gather(subtype, contribution, -c(signature, sumVar)) %>%
  separate(subtype, c("type", "motif"), sep="[.]") #%>%
  # filter(type!="C_T")
  # gather(signature, subtype)

sl1 <- merge(sl1, bg_rates, by=c("type", "motif")) %>%
  # mutate(fold=contribution/rate)
  mutate(pct=(contribution/rate-1)*100)
  # mutate(pct=ifelse(contribution > rate, (contribution/rate-1)*100, 1-rate/contribution))

sl1$Category2 <- plyr::mapvalues(sl1$type, orderedcats1a, orderedcats2a)
sl1$Category2 <- factor(sl1$Category2, levels=orderedcats2a)

sigloads_plot <- ggplot(sl1, aes(x=motif, y=pct, fill=signature))+
  geom_bar(stat="identity")+
  geom_hline(yintercept=0, linetype="dashed")+
  scale_y_continuous(breaks=seq(-100,300,50))+
  # scale_y_continuous(breaks=c(-log(seq(5,2)), 0, log(seq(2,5))), labels=c(seq(-5,-2), 0, seq(2,5)))+
  scale_fill_manual(values=iwhPalette5)+
  facet_grid(signature~Category2, scales="free_x")+
  ylab("% change in mutation rate")+
  theme_bw()+
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=16),
        strip.text.y=element_text(size=16),
        strip.text.x=element_text(size=16),
        axis.text.x=element_text(angle=90, vjust=0.5, size=14))
sigloads_plot
ggsave("/mnt/norbert/projects/ERV_mutation_hotspots/figs/pct_chg_sigloads.png", width=16, height=8)

sigloads_plot + theme(axis.title.y=element_blank(), axis.text.x=element_blank())
ggsave("/mnt/norbert/projects/ERV_mutation_hotspots/figs/pct_chg_sigloads_small.png", width=8, height=4)


# biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("org.Hs.eg.db")
library("GenomicFeatures")


txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
gene_ids <- as.list(org.Hs.egALIAS2EG)
mmr_genes <- gene_ids[c("MSH2", "MSH3", "MSH6", "MLH1", "MLH3", "PMS1", "PMS2")]

mmr_genes_gr <- transcriptsBy(txdb, by="gene")[unlist(mmr_genes)]
seqlengths(mmr_genes_gr) <- seqlengths(hg19IdeogramCyto)[names(seqlengths(mmr_genes_gr))]

nmf_W <- read.csv("/mnt/norbert/projects/ERV_mutation_hotspots/data/nmf_1Mb_W.csv")
names(nmf_W) <- c("ID", paste0("sig", 1:5))

nmf_W1 <- nmf_W %>%
  # dplyr::select(-sig3) %>%
  mutate(sumVar = rowSums(.[2:6])) %>%
  mutate_at(vars(starts_with("sig")), funs(./sumVar)) %>%
  separate(ID, c("chr", "start", "end")) %>%
  mutate(chr=paste0("chr", chr)) %>%
  mutate(start=as.numeric(start), end=as.numeric(end)) %>%
  gather(signature, contribution, sig1:sig5) %>%
  mutate(contribution=rollmedian(contribution, 5, fill=NA)) %>%
  arrange(chr, start) %>%
  filter(!is.na(contribution))

nmf_W1_thresh <- nmf_W1 %>% 
  group_by(signature) %>% 
  summarise(mean=mean(contribution), sd=sd(contribution)) %>% 
  mutate(threshold=mean+1.8*sd)

nmf_W1 <- merge(nmf_W1, nmf_W1_thresh, by="signature") %>% 
  mutate(peak=ifelse(contribution>threshold, TRUE, FALSE)) %>%
  mutate(size=ifelse(contribution>threshold, TRUE, FALSE))

nmf_W1_gr <- GRanges(nmf_W1)
seqlengths(nmf_W1_gr) <- seqlengths(hg19IdeogramCyto)[names(seqlengths(nmf_W1_gr))]

data(hg19IdeogramCyto, package = "biovizBase")
hg19 <- keepSeqlevels(hg19IdeogramCyto, paste0("chr", c(1:22)), pruning.mode="tidy")

# hg19_s1_chrs <- keepSeqlevels(hg19IdeogramCyto, paste0("chr", c(2, 3, 5, 7, 14)), pruning.mode="tidy")
# nmf_W1_gr_s1_chrs <- keepSeqlevels(nmf_W1_gr, paste0("chr", c(2, 3, 5, 7, 14)), pruning.mode="tidy")

# iwhPalette5 <- c("#9555b4", "#98be57", "#514143", "#c7624e", "#92b5b5")
# ggplot(hg19_s1_chrs) + 
#   layout_karyogram(cytoband = TRUE)+
#   layout_karyogram(nmf_W1_gr_s1_chrs, 
#                    geom = "point", 
#                    ylim = c(11, 31), 
#                    aes(x=start, y=contribution, colour=signature),
#                    size=0.8, alpha=0.8)+
#   scale_colour_manual(values = iwhPalette5)+
#   layout_karyogram(unlist(mmr_genes_gr), geom="rect", colour="red",
#                    # aes(x=start, y= -5, label=tx_id), 
#                    ylim = c(0,15))+
#   facet_grid(seqnames~., switch="y")+
#   ylab("signature contribution")+
#   # geom_alignment(facets = seqnames~.)+
#   theme_bw()+
#   theme(strip.text.y=element_text(size=16, angle=-90),
#         axis.text.y=element_blank())

p <- ggplot(hg19) + 
  layout_karyogram(cytoband = TRUE)+
  layout_karyogram(nmf_W1_gr, 
                   geom = "point", 
                   ylim = c(11, 51), 
                   aes(x=start, y=contribution, colour=signature, alpha=peak, size=size))+
  scale_colour_manual(values = iwhPalette5)+
  scale_size_manual(values = c(1,2))+
  scale_alpha_manual(values = c(0.3,1))+
  # scale_x_discrete(position = "top")+
  layout_karyogram(unlist(mmr_genes_gr), geom="rect", colour="red",
                   # aes(x=start, y= -5, label=tx_id), 
                   ylim = c(0,15))+
  facet_grid(as.numeric(gsub("chr", "", seqnames))~., switch="y")+
  ylab("signature contribution")+
  guides(fill=FALSE, alpha=FALSE)+
  # geom_alignment(facets = seqnames~.)+
  theme_classic()+
  theme(strip.text.y=element_text(size=16, angle=180),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank(),
        # axis.text.x.top = element_text(angle = 90, vjust=0.5, hjust=0),
        legend.position="top")

p
ggsave("/mnt/norbert/projects/ERV_mutation_hotspots/figs/chrs_all.png", width=7, height=17)

hg19.1.11 <- keepSeqlevels(hg19IdeogramCyto, paste0("chr", 1:11), pruning.mode="tidy")
nmf_W1_gr.1.11 <- keepSeqlevels(nmf_W1_gr, paste0("chr", 1:11), pruning.mode="tidy")
mmr_genes_gr.1.11 <- keepSeqlevels(mmr_genes_gr, paste0("chr", 1:11), pruning.mode="tidy")

p.sig1 <- ggplot(hg19.1.11) + 
  layout_karyogram(cytoband = TRUE)+
  layout_karyogram(nmf_W1_gr.1.11, 
                   geom = "point", 
                   ylim = c(11, 31), 
                   aes(x=start, y=contribution, colour=signature, alpha=peak),
                   size=3)+
  scale_colour_manual(values = iwhPalette5)+
  layout_karyogram(unlist(mmr_genes_gr.1.11), geom="rect", colour="red",
                   # aes(x=start, y= -5, label=tx_id), 
                   ylim = c(0,15))+
  facet_grid(as.numeric(gsub("chr", "", seqnames))~., switch="y")+
  ylab("signature contribution")+
  guides(fill=FALSE, alpha=FALSE)+
  # geom_alignment(facets = seqnames~.)+
  theme_classic()+
  theme(strip.text.y=element_text(size=16, angle=180),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position="bottom")

p.1.11
ggsave("/mnt/norbert/projects/ERV_mutation_hotspots/figs/chrs1.11.png", width=9.75, height=13)


p.1.11 <- ggplot(hg19.1.11) + 
  layout_karyogram(cytoband = TRUE)+
  layout_karyogram(nmf_W1_gr.1.11, 
                   geom = "point", 
                   ylim = c(11, 31), 
                   aes(x=start, y=contribution, colour=signature, alpha=peak),
                   size=3)+
  scale_colour_manual(values = iwhPalette5)+
  layout_karyogram(unlist(mmr_genes_gr.1.11), geom="rect", colour="red",
                   # aes(x=start, y= -5, label=tx_id), 
                   ylim = c(0,15))+
  facet_grid(as.numeric(gsub("chr", "", seqnames))~., switch="y")+
  ylab("signature contribution")+
  guides(fill=FALSE, alpha=FALSE)+
  # geom_alignment(facets = seqnames~.)+
  theme_classic()+
  theme(strip.text.y=element_text(size=16, angle=180),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position="bottom")

p.1.11
ggsave("/mnt/norbert/projects/ERV_mutation_hotspots/figs/chrs1.11.png", width=9.75, height=13)

hg19.12.22 <- keepSeqlevels(hg19IdeogramCyto, paste0("chr", 12:22), pruning.mode="tidy")
nmf_W1_gr.12.22 <- keepSeqlevels(nmf_W1_gr, paste0("chr", 12:22), pruning.mode="tidy")
mmr_genes_gr.12.22 <- keepSeqlevels(mmr_genes_gr, paste0("chr", 12:22), pruning.mode="tidy")


p.12.22 <- ggplot(hg19.12.22) + 
  layout_karyogram(cytoband = TRUE)+
  layout_karyogram(nmf_W1_gr.12.22, 
                   geom = "point", 
                   ylim = c(5, 33), 
                   aes(x=start, y=contribution, colour=signature, alpha=peak),
                   size=3)+
  scale_colour_manual(values = iwhPalette5)+
  layout_karyogram(unlist(mmr_genes_gr.12.22), geom="rect", colour="red",
                   # aes(x=start, y= -5, label=tx_id), 
                   ylim = c(0,15))+
  facet_grid(as.numeric(gsub("chr", "", seqnames))~., switch="y")+
  ylab("signature contribution")+
  guides(fill=FALSE, alpha=FALSE)+
  # geom_alignment(facets = seqnames~.)+
  theme_classic()+
  theme(strip.text.y=element_text(size=16, angle=180),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position="bottom")

p.12.22
ggsave("/mnt/norbert/projects/ERV_mutation_hotspots/figs/chrs12.22.png", width=9.75, height=13)

# clustchrs <- hg19
clustchrs<-hg19[seqnames(hg19) %in% unique(seqnames(clust_markers)),]
seqlengths(clustchrs) <- seqlengths(hg19Ideogram)[names(seqlengths(clustchrs))]

labs <- bins_nmf_long %>%
  group_by(CHR) %>%
  summarise(start=min(window),
    end=max(window),
    window=ceiling(start+(end-start)/2)) %>%
  mutate(fillcol=ifelse(CHR%%2==0, "even", "odd"))

iwhPalette5 <- c("#9555b4", "#98be57", "#514143", "#c7624e", "#92b5b5")
ggplot()+
  facet_wrap(~cluster, ncol=1, scales="free_y", strip.position="right")+
  geom_rect(data=labs, aes(xmin=start, xmax=end, ymin=0, ymax=1, fill=fillcol), alpha=0.5)+
  geom_hline(data=mean_contrib, aes(yintercept=mean_sig))+
  geom_hline(data=mean_contrib, aes(yintercept=mean_sig+2*sd_sig), linetype="dashed")+
  # geom_hline(data=mean_contrib, aes(yintercept=mean-2*sd), linetype="dashed")+
  geom_line(data=bins_nmf_long2,
    aes(x=window, y=prob_m, colour=cluster, group=cluster))+
  scale_fill_manual(values=c("grey70", "grey80"), guide=FALSE)+
  scale_colour_manual(values=iwhPalette5)+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(expand=c(0,0), breaks=labs$window, labels=labs$CHR)+
  xlab("chromosome")+
  ylab("signature contribution")+
  theme_bw()+
  theme(legend.position="bottom",
    axis.title.x=element_text(size=16),
    axis.title.y=element_text(size=16),
    strip.text.y=element_text(size=16),
    axis.text.x=element_text(size=14))
ggsave("figs/window_contributions.png", width=14, height=7)

# sigloads <- coef(bins_nmf) %>% 
#   data.frame %>% 
#   mutate(signature = paste0("X",row_number())) %>% 
#   mutate(sumVar = rowSums(.[1:96])) %>%
#   mutate_at(vars(contains("_")), funs(./sumVar)) %>%
#   gather(subtype, contribution, -c(signature, sumVar)) %>% 
#   separate(subtype, into=c("Type", "Motif"), 5) %>% 
#   mutate(Motif=substr(Motif,2,4))

sigloads$Type <- factor(sigloads$Type, levels=c("AT_CG", "AT_GC", "AT_TA", "GC_TA", "GC_CG", "GC_AT"))
levels(sigloads$Type) <- c("A_C", "A_G", "A_T", "C_A", "C_G", "C_T")

ggplot(sigloads, aes(x=Motif, y=contribution, fill=signature))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=iwhPalette5)+
  facet_grid(signature~Type, scales="free")+
  theme_bw()+
    theme(legend.position="bottom",
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=16),
    strip.text.y=element_text(size=16),
    strip.text.x=element_text(size=16),
    axis.text.x=element_text(angle=90, vjust=0.5, size=16))
ggsave("sigloads.png", width=14, height=7)  

sigloads_window <- sigloads %>% 
  mutate(subtype=paste0(Type, ".", Motif)) %>% 
  dplyr::select(signature, subtype, contribution) %>% 
  spread(signature, contribution)
sigloads_gene <- sig_loads %>% gather(subtype, contribution, -Sig) %>% spread(Sig, contribution)
sigloads_c <- merge(sigloads_window, sigloads_gene, by="subtype")
names(sigloads_c) <- c("Subtype", paste0("Window_Sig", 1:5), paste0("Gene_Sig", 1:5))

cs<-read.table("data/cancer_sigs.txt", header=T, sep="\t")
cs<-cs[,1:33]
snames<-paste0("COSMIC_Sig", 1:30)
names(cs)<-c("Type", "Seq3", "Subtype", snames)

cs$Type<-as.factor(cs$Type)
levels(cs$Type)<-c("C_A", "C_G", "C_T", "A_T", "A_G", "A_C")

revcomp <- function(DNAstr) {
  step1 <- chartr("ACGT","TGCA",DNAstr)
  step2 <- unlist(strsplit(step1, split=""))
  step3 <- rev(step2)
  step4 <- paste(step3, collapse="")
  return(step4)
}

cs$Seq3a<-ifelse(substr(cs$Type,1,1)=="A", as.character(reverse(complement(DNAStringSet(cs$Seq3)))),cs$Seq3)


cs2 <- cs %>%
  mutate(Motif = ifelse(substr(Type,1,1)=="A", 
                        as.character(reverse(complement(DNAStringSet(Seq3)))),
                        Seq3)) %>% 
  mutate(Subtype = paste0(Type, ".", Motif)) %>%
  dplyr::select(Subtype, COSMIC_Sig1:COSMIC_Sig30)

sigloads_c2 <- merge(sigloads_c, cs2, by="Subtype")

hm<-heatmapr(cor(sigloads_c[,2:11]), margins=c(100,100))
heatmaply(cor(sigloads_c[,2:11]), margins=c(100,100))


hm<-heatmapr(cor(sigloads_c2[,2:41]), margins=c(100,100))
heatmaply(cor(sigloads_c2[,2:41]), margins=c(100,100))

sigloads_c %>% gather(Signature, Contribution, -Subtype) %>% #head
  separate(Subtype, into=c("Type", "Motif"), sep="[.]") %>%
  mutate(Signature=factor(Signature, levels=rev(get_leaves_attr(hm$rows, "label")))) %>%
  ggplot(aes(x=Motif, y=Contribution, fill=Signature))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c(iwhPalette5, viridis(5))[c(10,3,8,5,9,6,2,1,4,7)])+
  facet_grid(Signature~Type, scales="free")+
  theme_bw()+
    theme(legend.position="bottom",
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=16),
    strip.text.y=element_text(size=16, angle=0),
    strip.text.x=element_text(size=16),
    axis.text.x=element_blank())
ggsave("figs/combined_sigloads.png", width=10, height=10)
