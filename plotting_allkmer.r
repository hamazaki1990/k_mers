library(tidyverse)
library(gridExtra)
library(grid)
setwd("~/github/k_mers/")

(peri<-read_csv("peri_and_livingcent.bed", col_names = c("chr", "start", "end")))
(cent<-read_csv("cent_annotation.csv")%>%filter(State=="live")%>%select(chr, CentStart, CentEnd))
(chromlen<-read_csv("chrom_length.csv")%>%select(chr, ChromStart, ChromEnd))
chromlen%>%left_join(cent)%>%mutate("qarm"=ChromEnd-CentEnd)%>%summarise(max(qarm))
chromlen%>%left_join(cent)%>%mutate("parm"=ChromStart-CentStart)%>%summarise(max(parm))

k = 18
hor = 11
(gfffile = str_c(k, "mer_chr", hor, "HOR_hit_test.gff",sep = ""))
(kmerhit<-kmerhit%>%mutate("kmermiddle"=(kmerstart+kmerend)/2))
(kmerhit <- read_csv(gfffile, col_names = c("seq", "source", "feature", "kmerstart", "kmerend", "score", "strand", "attribute"))
  %>% filter(feature == "kmer_hit") %>%separate(seq, c("filename", "chr"), ":")
  %>% select(chr, kmerstart, kmerend, strand))

(kmer_rel<-kmerhit%>%left_join(cent)%>%filter(kmermiddle<CentStart)%>%mutate("kmer_rel"=kmermiddle-CentStart)
  %>%bind_rows(kmerhit%>%left_join(cent)%>%filter(kmermiddle>CentEnd)%>%mutate("kmer_rel"=kmermiddle-CentEnd)))
(cnt_rel<-kmer_rel%>%filter(kmer_rel>0)%>%filter(strand=="+")%>%group_by(chr)%>%count(cut_width(kmer_rel, 100000, boundary=0))%>%ungroup%>%mutate(strand=1, strand_color="+")%>%
    bind_rows(kmer_rel%>%filter(kmer_rel>0)%>%filter(strand=="-")%>%group_by(chr)%>%count(cut_width(kmer_rel, 100000, boundary=0))%>%ungroup%>%mutate(strand=-1, strand_color="-"))%>%
    bind_rows(kmer_rel%>%filter(kmer_rel<0)%>%filter(strand=="+")%>%group_by(chr)%>%count(cut_width(kmer_rel, 100000, boundary=0))%>%ungroup%>%mutate(strand=1, strand_color="+"))%>%
    bind_rows(kmer_rel%>%filter(kmer_rel<0)%>%filter(strand=="-")%>%group_by(chr)%>%count(cut_width(kmer_rel, 100000, boundary=0))%>%ungroup%>%mutate(strand=-1, strand_color="-"))%>%
    mutate("position"=str_sub(`cut_width(kmer_rel, 1e+05, boundary = 0)`, 2, -2))%>%separate(position,c("kstart", "kend"),",")%>%
    mutate("count"=n*strand))

cnt_rel$kstart<-as.integer(cnt_rel$kstart)
cnt_rel$kend<-as.integer(cnt_rel$kend)
cnt_rel<-cnt_rel%>%mutate("x_rel"=(kstart+kend)/2)

cent_anno<-cent%>%select(chr)%>%mutate("annotation"="p_euchromatin", "plotboundary"=-5000000, "vlinealpha"=0)%>%
  bind_rows(cent%>%select(chr)%>%mutate("annotation"="p_pericentromere", "plotboundary"=0, "vlinealpha"=1))%>%
  bind_rows(cent%>%select(chr)%>%mutate("annotation"="q_pericentromere", "plotboundary"=0, "vlinealpha"=1))%>%
  bind_rows(cent%>%select(chr)%>%mutate("annotation"="q_euchromatin", "plotboundary"=5000000, "vlinealpha"=0))

df<-cnt_rel%>%filter(x_rel < -5000000)%>%mutate("annotation"="p_euchromatin")%>%
  bind_rows(cnt_rel%>%filter(x_rel<0 & x_rel>=-5000000)%>%mutate("annotation"="p_pericentromere"))%>%
  bind_rows(cnt_rel%>%filter(x_rel<5000000 & x_rel>=0)%>%mutate("annotation"="q_pericentromere"))%>%
  bind_rows(cnt_rel%>%filter(x_rel>=5000000)%>%mutate("annotation"="q_euchromatin"))%>%
  bind_rows(cent_anno)%>%
  mutate("cnt_log10"=strand*log10(n))
df
chromosomes<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX", "chrY")
df$chr_f = factor(df$chr, levels=chromosomes)
chrom_anno<-c("p_euchromatin", "p_pericentromere", "q_pericentromere", "q_euchromatin")
df$annotation_f = factor(df$annotation, levels=chrom_anno)
gp <-  ggplot(df) + geom_rect(aes(xmin=kstart, xmax=kend, ymin=0, ymax=count, fill=strand_color))
# gp <- ggplot(cnt_rel%>%filter(x_rel<5000000 & x_rel>-5000000)) + geom_col(aes(x=x_rel, y=count))
# gp<- gp + geom_rect(aes(xmin=CentStart, xmax=CentEnd, ymin=-Inf, ymax=Inf, alpha=0.5))
gp<- gp + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~annotation_f, scales="free_x", shrink = FALSE) # + geom_segment(aes(x=-5000000, xend=5000000, y=0, yend=0))
gp <- gp + scale_alpha_identity() + theme_bw()

print(gp+xlim(-10316944, -5000000))

gp_grobed1<-ggplotGrob(gp+xlim(-10316944, -5000000))
gp_grobed2<-ggplotGrob(gp+xlim(-5000000, 0))
gp_grobed3<-ggplotGrob(gp+xlim(0, 5000000))
gp_grobed4<-ggplotGrob(gp+theme(legend.position = "none")+xlim(5000000, 148102972))

gp_grobed$grobs[[3]]

gp_grobed$grobs$name

g_test<-gp_grobed

g_test$grobs[[102]] <-g_test$grobs[[103]]

grid.newpage()
grid.draw(g_test)
grid.newpage()
grid.draw(gp_grobed1)
grid.newpage()
grid.draw(gp_grobed2)
grid.newpage()
grid.draw(gp_grobed3)
grid.newpage()
grid.draw(gp_grobed4)

allkmer<-gp_grobed1

for (i in 26:49){
  allkmer$grobs[[i]]<-gp_grobed2$grobs[[i]]
}

for (i in 50:73){
  allkmer$grobs[[i]]<-gp_grobed3$grobs[[i]]
}

for (i in 74:97){
  allkmer$grobs[[i]]<-gp_grobed4$grobs[[i]]
}


allkmer$grobs[[103]]<-gp_grobed2$grobs[[103]]
allkmer$grobs[[104]]<-gp_grobed3$grobs[[104]]
allkmer$grobs[[105]]<-gp_grobed4$grobs[[105]]

grid.newpage()
grid.draw(allkmer)


cnt_rel%>%filter(kend<0)
df<-cnt_rel%>%filter(x_rel>=0)%>%mutate("kstart_log10"=log10(kstart), "kend_log10"=log10(kend), "xrel_log10"=log10(x_rel))%>%
  bind_rows(cnt_rel%>%filter(x_rel<0)%>%mutate("kstart_log10"=log10(abs(kstart)), "kend_log10"=log10(abs(kend)),"xrel_log10"=log10(abs(x_rel))))%>%
  mutate("cnt_log10"=strand*log10(n))
df%>%filter(x_rel<0)
chromosomes<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX", "chrY")
df$chr_f = factor(df$chr, levels=chromosomes)

gp1 <- ggplot(df%>%filter(x_rel>=0)) + geom_col(aes(x=xrel_log10, y=count, fill=strand_color))
gp1<- gp1 + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.) # + geom_segment(aes(x=-5000000, xend=5000000, y=0, yend=0))
gp1 <- gp1 + scale_alpha_identity() + theme_bw()
print(gp1)

gp2 <- ggplot(df%>%filter(x_rel<0)%>%mutate("kstart_ab"=-kstart, "kend_ab"=-kend)) + geom_rect(aes(xmin=kstart_ab, xmax=kend_ab, ymin=0, ymax=count, fill=strand_color))
gp2<- gp2 + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.) # + geom_segment(aes(x=-5000000, xend=5000000, y=0, yend=0))
gp2 <- gp2 + scale_alpha_identity() + scale_x_reverse() + scale_x_log10() + theme_bw()
print(gp2)




k = 18
hor = 11
(gfffile = str_c(k, "mer_chr", hor, "HOR_hit_test.gff",sep = ""))
(kmerhit<-kmerhit%>%mutate("kmermiddle"=(kmerstart+kmerend)/2))
(kmerhit <- read_csv(gfffile, col_names = c("seq", "source", "feature", "kmerstart", "kmerend", "score", "strand", "attribute"))
  %>% filter(feature == "kmer_hit") %>%separate(seq, c("filename", "chr"), ":")
  %>% select(chr, kmerstart, kmerend, strand))

(kmer_rel<-kmerhit%>%left_join(cent)%>%filter(kmermiddle<CentStart)%>%mutate("kmer_rel"=kmermiddle-CentStart)
  %>%bind_rows(kmerhit%>%left_join(cent)%>%filter(kmermiddle>CentEnd)%>%mutate("kmer_rel"=kmermiddle-CentEnd)))
(cnt_rel<-kmer_rel%>%filter(kmer_rel>0)%>%filter(strand=="+")%>%group_by(chr)%>%count(cut_width(kmer_rel, 100000, boundary=0))%>%ungroup%>%mutate(strand=1, strand_color="+")%>%
    bind_rows(kmer_rel%>%filter(kmer_rel>0)%>%filter(strand=="-")%>%group_by(chr)%>%count(cut_width(kmer_rel, 100000, boundary=0))%>%ungroup%>%mutate(strand=-1, strand_color="-"))%>%
    bind_rows(kmer_rel%>%filter(kmer_rel<0)%>%filter(strand=="+")%>%group_by(chr)%>%count(cut_width(kmer_rel, 100000, boundary=0))%>%ungroup%>%mutate(strand=1, strand_color="+"))%>%
    bind_rows(kmer_rel%>%filter(kmer_rel<0)%>%filter(strand=="-")%>%group_by(chr)%>%count(cut_width(kmer_rel, 100000, boundary=0))%>%ungroup%>%mutate(strand=-1, strand_color="-"))%>%
    mutate("position"=str_sub(`cut_width(kmer_rel, 1e+05, boundary = 0)`, 2, -2))%>%separate(position,c("kstart", "kend"),",")%>%
    mutate("count"=n*strand))

cnt_rel$kstart<-as.integer(cnt_rel$kstart)
cnt_rel$kend<-as.integer(cnt_rel$kend)
cnt_rel<-cnt_rel%>%mutate("x_rel"=(kstart+kend)/2)

cent_anno<-cent%>%select(chr)%>%mutate("annotation"="p_euchromatin", "plotboundary"=-5000000, "vlinealpha"=0)%>%
  bind_rows(cent%>%select(chr)%>%mutate("annotation"="p_pericentromere", "plotboundary"=0, "vlinealpha"=1))%>%
  bind_rows(cent%>%select(chr)%>%mutate("annotation"="q_pericentromere", "plotboundary"=0, "vlinealpha"=1))%>%
  bind_rows(cent%>%select(chr)%>%mutate("annotation"="q_euchromatin", "plotboundary"=5000000, "vlinealpha"=0))

df<-cnt_rel%>%filter(x_rel < -5000000)%>%mutate("annotation"="p_euchromatin")%>%
  bind_rows(cnt_rel%>%filter(x_rel<0 & x_rel>=-5000000)%>%mutate("annotation"="p_pericentromere"))%>%
  bind_rows(cnt_rel%>%filter(x_rel<5000000 & x_rel>=0)%>%mutate("annotation"="q_pericentromere"))%>%
  bind_rows(cnt_rel%>%filter(x_rel>=5000000)%>%mutate("annotation"="q_euchromatin"))%>%
  bind_rows(cent_anno)%>%
  mutate("cnt_log10"=strand*log10(n))
df
chromosomes<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX", "chrY")
df$chr_f = factor(df$chr, levels=chromosomes)
chrom_anno<-c("p_euchromatin", "p_pericentromere", "q_pericentromere", "q_euchromatin")
df$annotation_f = factor(df$annotation, levels=chrom_anno)
gp <-  ggplot(df) + geom_rect(aes(xmin=kstart, xmax=kend, ymin=0, ymax=count, fill=strand_color))
# gp <- ggplot(cnt_rel%>%filter(x_rel<5000000 & x_rel>-5000000)) + geom_col(aes(x=x_rel, y=count))
# gp<- gp + geom_rect(aes(xmin=CentStart, xmax=CentEnd, ymin=-Inf, ymax=Inf, alpha=0.5))
gp<- gp + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.) # + geom_segment(aes(x=-5000000, xend=5000000, y=0, yend=0))
gp <- gp + scale_alpha_identity() + theme_bw() 

gp1 <- ggplot(df%>%filter(annotation=="p_euchromatin")) + geom_rect(aes(xmin=kstart, xmax=kend, ymin=0, ymax=count, fill=strand_color))
gp1 <-gp1 + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.) + scale_alpha_identity() + theme_bw() 
gp1 <-gp1 + coord_cartesian(xlim=c(-10316944, -5000000)) + theme(strip.background = element_blank(), strip.text.y = element_blank(),legend.position = "none")

gp2 <- ggplot(df%>%filter(annotation=="p_pericentromere")) + geom_rect(aes(xmin=kstart, xmax=kend, ymin=0, ymax=count, fill=strand_color))
gp2 <-gp2 + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.) + scale_alpha_identity() + theme_bw() 
gp2 <-gp2 + coord_cartesian(xlim=c(-5000000, 0)) + theme(strip.background = element_blank(), strip.text.y = element_blank(),legend.position = "none")

gp3 <- ggplot(df%>%filter(annotation=="q_pericentromere")) + geom_rect(aes(xmin=kstart, xmax=kend, ymin=0, ymax=count, fill=strand_color))
gp3 <-gp3 + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.) + scale_alpha_identity() + theme_bw() 
gp3 <-gp3 + coord_cartesian(xlim=c(0, 5000000)) + theme(strip.background = element_blank(), strip.text.y = element_blank(),legend.position = "none")


print(gp2)





gp_grobed1<-ggplotGrob(gp+xlim(-10316944, -5000000))
gp_grobed2<-ggplotGrob(gp+xlim(-5000000, 0))
gp_grobed3<-ggplotGrob(gp+xlim(0, 5000000))
gp_grobed4<-ggplotGrob(gp+theme(legend.position = "none")+xlim(5000000, 148102972))

















ggsave("chr11_18mer_peri.png", gp, width=10, height=10)