library(tidyverse)
setwd("~/github/k_mers/")

(peri<-read_csv("peri_and_livingcent.bed", col_names = c("chr", "start", "end")))
(cent<-read_csv("cent_annotation.csv")%>%filter(State=="live")%>%select(chr, CentStart, CentEnd))

k = 18
hor = 11
(gfffile = str_c(k, "mer_chr", hor, "HOR_hit_test.gff",sep = ""))
(kmerhit <- read_csv(gfffile, col_names = c("seq", "source", "feature", "kmerstart", "kmerend", "score", "strand", "attribute"))
          %>% filter(feature == "kmer_hit") %>%separate(seq, c("filename", "chr"), ":")
          %>% select(chr, kmerstart, kmerend, strand))
nclass.Sturges(kmerhit)

(kmerhit<-kmerhit%>%mutate("kmermiddle"=(kmerstart+kmerend)/2))
(cnt<-kmerhit%>%filter(strand=="+")%>%group_by(chr)%>%count(cut_width(kmermiddle, 100000))%>%ungroup%>%mutate(strand=1)
  %>%bind_rows(kmerhit%>%filter(strand=="-")%>%group_by(chr)%>%count(cut_width(kmermiddle, 100000))%>%ungroup%>%mutate(strand=-1))
  %>%mutate("position"=str_sub(`cut_width(kmermiddle, 1e+05)`, 2, -2))%>%separate(position,c("kstart", "kend"),",")
  %>%mutate("count"=n*strand))

cnt$kstart<-as.integer(cnt$kstart)
cnt$kend<-as.integer(cnt$kend)
(cnt<-cnt%>%mutate("x"=(kstart+kend)/2))
# cnt%>%View
(df<-cnt%>%bind_rows(cent)%>%filter(chr=="chr1"))
df
gp <- ggplot(df) + geom_col(aes(x=x, y=count))
gp <- gp + geom_rect(aes(xmin=kstart, xmax=kend, ymin=min(count, 0), ymax=max(count,0)))
gp <- gp + geom_rect(aes(xmin=CentStart, xmax=CentEnd, ymin=-Inf, ymax=Inf, alpha=0.5))
gp <- gp + xlim(cent$CentStart[1]-5000000,cent$CentEnd[1]+5000000)
# gp<- gp + facet_grid(chr~.,scales="free")
print(gp)


(df<-bind_rows(cnt,cent))
chromosomes<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX", "chrY")
df$chr_f = factor(df$chr, levels=chromosomes)
gp <- ggplot(df) + geom_col(aes(x=x, y=count))
gp<- gp + geom_rect(aes(xmin=CentStart, xmax=CentEnd, ymin=-Inf, ymax=Inf, alpha=0.5))
df
gp <- gp + facet_wrap(~chr_f,scales="free")
gp <- gp + theme_minimal()
print(gp)


plotlim<-cent%>%mutate(PeriStart=CentStart-3000000, PeriEnd=CentEnd+3000000)%>%select(chr, PeriStart, PeriEnd)%>%gather(PeriStart,PeriEnd,key="plotlim", value="x")
plotlim
cnt
(df<-bind_rows(cnt%>%mutate("alpha"=1),plotlim%>%select(chr, x)%>%mutate("alpha"=0, "count"=0),cent))
# df%>%View
chromosomes<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX", "chrY")
df$chr_f = factor(df$chr, levels=chromosomes)
gp <- ggplot(df) + geom_col(aes(x=x, y=count, alpha=alpha))
gp<- gp + geom_rect(aes(xmin=CentStart, xmax=CentEnd, ymin=-Inf, ymax=Inf, alpha=0.5))
gp<- gp + facet_wrap(~chr_f,scales="free")
print(gp)

(kmer_rel<-kmerhit%>%left_join(cent)%>%filter(kmermiddle<CentStart)%>%mutate("kmer_rel"=kmermiddle-CentStart)
  %>%bind_rows(kmerhit%>%left_join(cent)%>%filter(kmermiddle>CentEnd)%>%mutate("kmer_rel"=kmermiddle-CentEnd)))
(cnt_rel<-kmer_rel%>%filter(kmer_rel>0)%>%filter(strand=="+")%>%group_by(chr)%>%count(cut_width(kmer_rel, 100000, boundary=0))%>%ungroup%>%mutate(strand=1, strand_color="+")%>%
  bind_rows(kmer_rel%>%filter(kmer_rel>0)%>%filter(strand=="-")%>%group_by(chr)%>%count(cut_width(kmer_rel, 100000, boundary=0))%>%ungroup%>%mutate(strand=-1, strand_color="-"))%>%
  bind_rows(kmer_rel%>%filter(kmer_rel<0)%>%filter(strand=="+")%>%group_by(chr)%>%count(cut_width(kmer_rel, 100000, boundary=0))%>%ungroup%>%mutate(strand=1, strand_color="+"))%>%
  bind_rows(kmer_rel%>%filter(kmer_rel<0)%>%filter(strand=="-")%>%group_by(chr)%>%count(cut_width(kmer_rel, 100000, boundary=0))%>%ungroup%>%mutate(strand=-1, strand_color="-"))%>%
  mutate("position"=str_sub(`cut_width(kmer_rel, 1e+05, boundary = 0)`, 2, -2))%>%separate(position,c("kstart", "kend"),",")%>%
  mutate("count"=n*strand))

#kmer_rel%>%filter(kmer_rel<100000 &kmer_rel>-100000)%>%View
cnt_rel$kstart<-as.integer(cnt_rel$kstart)
cnt_rel$kend<-as.integer(cnt_rel$kend)
cnt_rel<-cnt_rel%>%mutate("x_rel"=(kstart+kend)/2)

# df<-cnt_rel%>%bind_rows(cent)%>%filter(chr=="chr1")
# df
# gp <- ggplot(df) + geom_col(aes(x=x_rel, y=count))
# print(gp)

cnt_rel

cent_anno<-cent%>%select(chr)%>%mutate("annotation"="p_euchromatin", "plotboundary"=-5000000, "vlinealpha"=0)%>%
  bind_rows(cent%>%select(chr)%>%mutate("annotation"="p_pericentromere", "plotboundary"=0, "vlinealpha"=1))%>%
  bind_rows(cent%>%select(chr)%>%mutate("annotation"="q_pericentromere", "plotboundary"=0, "vlinealpha"=1))%>%
  bind_rows(cent%>%select(chr)%>%mutate("annotation"="q_euchromatin", "plotboundary"=5000000, "vlinealpha"=0))

df<-cnt_rel%>%filter(x_rel < -5000000)%>%mutate("annotation"="p_euchromatin")%>%
  bind_rows(cnt_rel%>%filter(x_rel<0 & x_rel>=-5000000)%>%mutate("annotation"="p_pericentromere"))%>%
  bind_rows(cnt_rel%>%filter(x_rel<5000000 & x_rel>=0)%>%mutate("annotation"="q_pericentromere"))%>%
  bind_rows(cnt_rel%>%filter(x_rel>=5000000)%>%mutate("annotation"="q_euchromatin"))%>%
  bind_rows(cent_anno)

chromosomes<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX", "chrY")
df$chr_f = factor(df$chr, levels=chromosomes)
chrom_anno<-c("p_euchromatin", "p_pericentromere", "q_pericentromere", "q_euchromatin")
df$annotation_f = factor(df$annotation, levels=chrom_anno)
gp <- ggplot(df) + geom_rect(aes(xmin=kstart, xmax=kend, ymin=0, ymax=count, fill=strand_color))
# gp <- ggplot(cnt_rel%>%filter(x_rel<5000000 & x_rel>-5000000)) + geom_col(aes(x=x_rel, y=count))
# gp<- gp + geom_rect(aes(xmin=CentStart, xmax=CentEnd, ymin=-Inf, ymax=Inf, alpha=0.5))
gp<- gp + geom_segment(aes(x=plotboundary, xend=plotboundary, y=-Inf, yend=Inf, alpha=vlinealpha), colour="#999999") + facet_grid(chr_f~annotation_f,scales="free") # + geom_segment(aes(x=-5000000, xend=5000000, y=0, yend=0))
gp <- gp + scale_alpha_identity() + theme_bw()
print(gp)
ggsave("chr11_18mer_peri.png", gp, width=10, height=10)

peri

cnt%>%filter(x>130000000)


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
hor = 11
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

k = 19
hor = 11
(gfffile = str_c(k, "mer_chr", hor, "HOR_hit_test.gff",sep = ""))
(kmeratt <- read_csv(gfffile, col_names = c("seq", "source", "feature", "kmerstart", "kmerend", "score", "strand", "attribute"))
  %>% filter(feature == "kmer_hit") %>%separate(seq, c("filename", "chr"), ":")
  %>% select(chr, kmerstart, kmerend, strand, attribute))

kmeratt%>%filter(strand=="+")%>%group_by(attribute)%>%mutate("n"=n())%>%arrange(-n)%>%ungroup()%>%select(attribute, n)%>%distinct()%>%View
kmeratt%>%filter(strand=="-")%>%group_by(attribute)%>%mutate("n"=n())%>%arrange(-n)%>%ungroup()%>%select(attribute, n)%>%distinct()%>%View

k = 18
hor = 1
(kmerfile = str_c("chr", hor, "_", k, "merlist.csv",sep = ""))
(kmerlist <- read_csv(kmerfile, col_names = c("index", "seq")))
for (hor in c(2, 11)){
  (kmerfile = str_c("chr", hor, "_", k, "merlist.csv",sep = ""))
  (l<-read_csv(kmerfile, col_names = c("index", "seq")))
  (kmerlist <- kmerlist%>%full_join(l, by="seq"))
}
kmerlist%>%filter(is.na(index.x)==FALSE & is.na(index)==FALSE|is.na(index.y)==FALSE & is.na(index)==FALSE|is.na(index.x)==FALSE & is.na(index.y)==FALSE)%>%View
