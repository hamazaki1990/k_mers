library(tidyverse)
getwd()
setwd("github/k_mers/")

(ch<-read_csv("chrom_length.csv"))
(ch<-ch%>%mutate("10mer"=ChromEnd/(4**10), "12mer"=ChromEnd/(4**12), "14mer"=ChromEnd/(4**14), "15mer"=ChromEnd/(4**15), "16mer"=ChromEnd/(4**16)))

(hit<-read.csv("chr2_20mer_chr1HOR_hit_test.csv", col_names = c("cent", "kmer", "seq", "chr", "position")))

gc <- 