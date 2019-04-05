library(tidyverse)
setwd("~/github/k_mers/")

(peri<-read_csv("peri_and_livingcent.bed", col_names = c("chr", "start", "end")))
(cent<-read_csv("cent_annotation.csv")%>%filter(State=="live")%>%select(chr, CentStart, CentEnd))

k = 16
hor = 1
(gfffile = str_c(k, "mer_chr", hor, "HOR_hit_test.gff",sep = ""))
(kmerhit <- read_csv(gfffile, col_names = c("seq", "source", "feature", "kmerstart", "kmerend", "score", "strand", "attribute"))
          %>% filter(feature == "kmer_hit") %>%separate(seq, c("filename", "chr"), ":")
          %>% select(chr, kmerstart, kmerend, strand))

(cnt<-kmerhit%>%left_join(peri, by="chr")%>%filter(start<=kmerstart & kmerend<=end)
  %>%filter(strand=="+")%>%group_by(chr)%>%count(cut_width(kmerstart, 1000))%>%ungroup%>%mutate(strand=1)
  %>%bind_rows(kmerhit%>%filter(strand=="-")%>%group_by(chr)%>%count(cut_width(kmerstart, 1000))%>%ungroup%>%mutate(strand=-1))
  %>%mutate("position"=str_sub(`cut_width(kmerstart, 1000)`, 2, -2))%>%separate(position,c("kstart", "kend"),",")
  %>%mutate("count"=n*strand))

cnt$kstart<-as.integer(cnt$kstart)
cnt$kend<-as.integer(cnt$kend)
(cnt<-cnt%>%mutate("x"=(kstart+kend)/2))

(df<-cnt%>%bind_rows(cent)%>%filter(chr=="chr1"))
gp <- ggplot(df) + geom_col(aes(x=x, y=count))
gp<- gp + geom_rect(aes(xmin=CentStart, xmax=CentEnd, ymin=-Inf, ymax=Inf, alpha=0.5))
gp<- gp + xlim(cent$CentStart[1],129785432)
# gp<- gp + facet_grid(chr~.,scales="free")
print(gp)

(df<-cnt%>%bind_rows(cent)%>%filter(chr=="chr2"))
gp <- ggplot(df) + geom_col(aes(x=x, y=count))
gp<- gp + geom_rect(aes(xmin=CentStart, xmax=CentEnd, ymin=-Inf, ymax=Inf, alpha=0.5))
gp<- gp + xlim(87188144,99090557)
# gp<- gp + facet_grid(chr~.,scales="free")
print(gp)

(df<-cnt%>%bind_rows(cent)%>%filter(chr=="chr3"))
gp <- ggplot(df) + geom_col(aes(x=x, y=count))
gp<- gp + geom_rect(aes(xmin=CentStart, xmax=CentEnd, ymin=-Inf, ymax=Inf, alpha=0.5))
gp<- gp + xlim(86553418,98655574)
# gp<- gp + facet_grid(chr~.,scales="free")
print(gp)

3000000000*(0.2^0.4)*(0.3^0.6)

hor = 1
k = 16
listfile = str_c("chr", hor, "_", k, "merlist.csv",sep = "")
sequence = str_c("HOR", hor, "_", k, "mer",sep = "")
(l <-read_csv(listfile, col_names = c("attribute", sequence)))
(gfffile = str_c(k, "mer_chr", hor, "HOR_hit_test.gff",sep = ""))
(kmeratt <- read_csv(gfffile, col_names = c("seq", "source", "feature", "kmerstart", "kmerend", "score", "strand", "attribute"))
  %>% filter(feature == "kmer_hit") %>%separate(seq, c("filename", "chr"), ":")
  %>% select(chr, kmerstart, kmerend, strand, attribute))
(HOR1_plus<-kmeratt%>%filter(strand=="+")%>%group_by(attribute)%>%mutate("n"=n())%>%arrange(-n)%>%ungroup()%>%select(attribute, n)%>%distinct()%>%
  left_join(l, by="attribute"))
for (k in 17:20){
  listfile = str_c("chr", hor, "_", k, "merlist.csv",sep = "")
  sequence = str_c("HOR", hor, "_", k, "mer",sep = "")
  (l <-read_csv(listfile, col_names = c("attribute", sequence)))
  (gfffile = str_c(k, "mer_chr", hor, "HOR_hit_test.gff",sep = ""))
  (kmeratt <- read_csv(gfffile, col_names = c("seq", "source", "feature", "kmerstart", "kmerend", "score", "strand", "attribute"))
    %>% filter(feature == "kmer_hit") %>%separate(seq, c("filename", "chr"), ":")
    %>% select(chr, kmerstart, kmerend, strand, attribute))
  (kmerplus<-kmeratt%>%filter(strand=="+")%>%group_by(attribute)%>%mutate("n"=n())%>%arrange(-n)%>%ungroup()%>%select(attribute, n)%>%distinct()
    %>%left_join(l, by="attribute"))
  HOR1_plus<-HOR1_plus%>%full_join(kmerplus, by="attribute")
}
HOR1_plus<-HOR1_plus%>%separate(attribute, c("HOR", "index"), ":")
HOR1_plus$index<-as.integer(HOR1_plus$index)
(HOR1_plus%>%arrange(index)%>%mutate("n.x_sum"=n.x+n.x%>%lead())%>%mutate(HOR1_16_17=n.y<n.x_sum)%>%filter(HOR1_16_17==FALSE))     # %>%write_csv("HOR1_kmerhit.csv", col_names = TRUE)
(HOR1_plus%>%arrange(index)%>%mutate("n.y_sum"=n.y+n.y%>%lead())%>%mutate(HOR1_17_18=n.x.x<n.y_sum)%>%filter(HOR1_17_18==FALSE))
(HOR1_plus%>%arrange(index)%>%mutate("n.x.x_sum"=n.x.x+n.x.x%>%lead())%>%mutate(HOR1_18_19=n.y.y<n.x.x_sum)%>%filter(HOR1_18_19==FALSE))

kmeratt%>%filter(strand=="-")%>%group_by(attribute)%>%mutate("n"=n())%>%arrange(-n)%>%ungroup()%>%select(attribute, n)%>%distinct()%>%View

k = 19
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
chr_kmerhit%>%View
chr_kmerhit%>%filter(endcheck==TRUE | startcheck==TRUE)%>%
  mutate("gap"==)
  View()

chr_kmerhit%>%filter(endcheck==TRUE | startcheck==TRUE)%>%mutate("kmergroup"=kmerstart-index)%>%View



kmeratt%>%filter(chr=="chr1")%>%filter(strand=="+")%>%group_by(attribute)%>%mutate("n"=n())%>%ungroup()%>%select(attribute, n)%>%distinct()%>%arrange(-n)%>%View
kmeratt%>%filter(strand=="-")%>%group_by(attribute)%>%mutate("n"=n())%>%arrange(-n)%>%ungroup()%>%select(attribute, n)%>%distinct()%>%View














k = 19
hor = 11
(gfffile = str_c(k, "mer_chr", hor, "HOR_hit_test.gff",sep = ""))
(kmeratt <- read_csv(gfffile, col_names = c("seq", "source", "feature", "kmerstart", "kmerend", "score", "strand", "attribute"))
  %>% filter(feature == "kmer_hit") %>%separate(seq, c("filename", "chr"), ":")
  %>% select(chr, kmerstart, kmerend, strand, attribute))

kmeratt%>%filter(strand=="+")%>%group_by(attribute)%>%mutate("n"=n())%>%arrange(-n)%>%ungroup()%>%select(attribute, n)%>%distinct()%>%View
kmeratt%>%filter(strand=="-")%>%group_by(attribute)%>%mutate("n"=n())%>%arrange(-n)%>%ungroup()%>%select(attribute, n)%>%distinct()%>%View

