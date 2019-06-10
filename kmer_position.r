library(tidyverse)
library(circlize)
install.packages("circlize")
setwd("~/github/k_mers/")

(peri<-read_csv("peri_and_livingcent.bed", col_names = c("chr", "start", "end")))
(cent<-read_csv("cent_annotation.csv")%>%filter(State=="live")%>%select(chr, CentStart, CentEnd))
(hg38gtf<-read_tsv("Homo_sapiens.GRCh38.90.gtf", col_names = c("chr", "sourse", "feature", "gtfstart", "gtfend", "score", "strand", "frame", "gtfattribute")))



k = 19
hor = 1
(gfffile = str_c("HOR", hor, "_", k, "mer_chr", hor, "HOR_hit_test.gff",sep = ""))
(kmeratt <- read_csv(gfffile, col_names = c("chr", "source", "feature", "kmerstart", "kmerend", "score", "strand", "attribute"), skip = 1)
  %>% filter(feature == "kmer_hit") %>%separate(attribute, c("HOR", "index"), ":")
  %>% select(chr, kmerstart, kmerend, strand, HOR, index))

kmeratt$index<-as.integer(kmeratt$index)

#gp<-ggplot(kmeratt)+geom_bar(aes(x=index))+scale_x_continuous(breaks = c(169, 340, 508))
#print(gp)

hor = 2
(gfffile = str_c("HOR", hor, "_", k, "mer_chr", hor, "HOR_hit_test.gff",sep = ""))
(kmeratt <- read_csv(gfffile, col_names = c("chr", "source", "feature", "kmerstart", "kmerend", "score", "strand", "attribute"), skip = 1)
  %>% filter(feature == "kmer_hit") %>%separate(attribute, c("HOR", "index"), ":")
  %>% select(chr, kmerstart, kmerend, strand, HOR, index))

kmeratt$index<-as.integer(kmeratt$index)

#gp2<-ggplot(kmeratt)+geom_bar(aes(x=index))+scale_x_continuous(breaks = c(171, 342, 513, 680, 851, 1022, 1193, 1361))
#print(gp2)


hor = 11
(gfffile = str_c("HOR", hor, "_", k, "mer_chr", hor, "HOR_hit_test.gff",sep = ""))
(kmeratt <- read_csv(gfffile, col_names = c("chr", "source", "feature", "kmerstart", "kmerend", "score", "strand", "attribute"), skip = 1)
  %>% filter(feature == "kmer_hit") %>%separate(attribute, c("HOR", "index"), ":")
  %>% select(chr, kmerstart, kmerend, strand, HOR, index))

kmeratt$index<-as.integer(kmeratt$index)

#gp3<-ggplot(kmeratt)+geom_bar(aes(x=index))+scale_x_continuous(breaks = c(170, 341, 509, 680))
#print(gp3)

#g<-plot_grid(gp, gp2, gp3, ncol=3)
#print(g)

#ggsave("HOR_HOR_19mer.png", g, width=10, height=6)
