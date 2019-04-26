library(tidyverse)
setwd("~/github/k_mers/")

(peri<-read_csv("peri_and_livingcent.bed", col_names = c("chr", "start", "end")))
(cent<-read_csv("cent_annotation.csv")%>%filter(State=="live")%>%select(chr, CentStart, CentEnd))
(hg38gtf<-read_tsv("Homo_sapiens.GRCh38.90.gtf", col_names = c("chr", "sourse", "feature", "gtfstart", "gtfend")))


k = 14
hor = 1
(gfffile = str_c(k, "mer_chr", hor, "HOR_hit_test.gff",sep = ""))
(kmeratt <- read_csv(gfffile, col_names = c("seq", "source", "feature", "kmerstart", "kmerend", "score", "strand", "attribute"))
  %>% filter(feature == "kmer_hit") %>%separate(seq, c("filename", "chr"), ":")
  %>% select(chr, kmerstart, kmerend, strand, attribute))
chr_kmerhit<-kmeratt%>%filter(chr=="chr1", strand=="+")%>%arrange(kmerstart)%>%
  #  mutate("kmerend_lead"=kmerend%>%lead())%>%mutate("kmerend_lag"=kmerend%>%lag())%>%
  mutate("startcheck"=((kmerstart-kmerstart%>%lag())!=1))%>%
  mutate("endcheck"=((kmerend%>%lead()-kmerend)!=1))%>%
  separate(attribute, c("HOR", "index"), ":")

chr_kmerhit$index<-as.integer(chr_kmerhit$index)
chr_kmerhit[is.na(chr_kmerhit)]<-TRUE
chr_kmerhit
chr_kmerhit%>%filter(endcheck==TRUE | startcheck==TRUE)%>%
  mutate("gap"=kmerstart%>%lead-kmerend)%>%
  #  filter(endcheck==TRUE)%>%
  mutate("kmergroup"=kmerstart-index)%>%
  View()

chr_kmerhit%>%filter(endcheck==TRUE | startcheck==TRUE)%>%mutate("kmergroup"=kmerstart-index)%>%View

kmeratt%>%filter(chr=="chr1")%>%filter(strand=="+")%>%group_by(attribute)%>%mutate("n"=n())%>%ungroup()%>%select(attribute, n)%>%distinct()%>%arrange(-n)%>%View
kmeratt%>%filter(strand=="-")%>%group_by(attribute)%>%mutate("n"=n())%>%arrange(-n)%>%ungroup()%>%select(attribute, n)%>%distinct()%>%View
