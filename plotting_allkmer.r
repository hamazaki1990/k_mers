library(tidyverse)
library(gridExtra)
library(grid)
library(cowplot)
setwd("~/github/k_mers/")

(peri<-read_csv("peri_and_livingcent.bed", col_names = c("chr", "start", "end")))
(cent<-read_csv("cent_annotation.csv")%>%filter(State=="live")%>%select(chr, CentStart, CentEnd))
(chromlen<-read_csv("chrom_length.csv")%>%select(chr, ChromStart, ChromEnd))
chromlen%>%left_join(cent)%>%mutate("qarm"=ChromEnd-CentEnd)%>%summarise(max(qarm))
chromlen%>%left_join(cent)%>%mutate("parm"=ChromStart-CentStart)%>%summarise(max(parm))

k = 18
hor = 11
(gfffile = str_c(k, "mer_chr", hor, "HOR_hit_test.gff",sep = ""))

(kmerhit <- read_csv(gfffile, col_names = c("seq", "source", "feature", "kmerstart", "kmerend", "score", "strand", "attribute"))
  %>% filter(feature == "kmer_hit") %>%separate(seq, c("filename", "chr"), ":")
  %>% select(chr, kmerstart, kmerend, strand))
(kmerhit<-kmerhit%>%mutate("kmermiddle"=(kmerstart+kmerend)/2))

(kmer_rel<-kmerhit%>%left_join(cent)%>%filter(kmermiddle<CentStart)%>%mutate("kmer_rel"=kmermiddle-CentStart)
  %>%bind_rows(kmerhit%>%left_join(cent)%>%filter(kmermiddle>CentEnd)%>%mutate("kmer_rel"=kmermiddle-CentEnd)))
kmer_rel<-kmer_rel%>%group_by(chr)%>%distinct(kmermiddle, .keep_all=TRUE)%>%ungroup
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

#gp_grobed$grobs[[3]]

#gp_grobed$grobs$name

#g_test<-gp_grobed

#g_test$grobs[[102]] <-g_test$grobs[[103]]

#grid.newpage()
#grid.draw(g_test)
#grid.newpage()
#grid.draw(gp_grobed1)
#grid.newpage()
#grid.draw(gp_grobed2)
#grid.newpage()
#grid.draw(gp_grobed3)
#grid.newpage()
#grid.draw(gp_grobed4)

#allkmer<-gp_grobed1

#for (i in 26:49){
#  allkmer$grobs[[i]]<-gp_grobed2$grobs[[i]]
#}

#for (i in 50:73){
#  allkmer$grobs[[i]]<-gp_grobed3$grobs[[i]]
#}

#for (i in 74:97){
#  allkmer$grobs[[i]]<-gp_grobed4$grobs[[i]]
#}


#allkmer$grobs[[103]]<-gp_grobed2$grobs[[103]]
#allkmer$grobs[[104]]<-gp_grobed3$grobs[[104]]
#allkmer$grobs[[105]]<-gp_grobed4$grobs[[105]]

#grid.newpage()
#grid.draw(allkmer)


k = 18
hor = 11
(gfffile = str_c(k, "mer_chr", hor, "HOR_hit_test.gff",sep = ""))
(kmerhit <- read_csv(gfffile, col_names = c("seq", "source", "feature", "kmerstart", "kmerend", "score", "strand", "attribute"))
  %>% filter(feature == "kmer_hit") %>%separate(seq, c("filename", "chr"), ":")
  %>% select(chr, kmerstart, kmerend, strand))
(kmerhit<-kmerhit%>%mutate("kmermiddle"=(kmerstart+kmerend)/2))


(kmer_rel<-kmerhit%>%left_join(cent)%>%filter(kmermiddle<CentStart)%>%mutate("kmer_rel"=kmermiddle-CentStart)
  %>%bind_rows(kmerhit%>%left_join(cent)%>%filter(kmermiddle>CentEnd)%>%mutate("kmer_rel"=kmermiddle-CentEnd)))
kmer_rel<-kmer_rel%>%group_by(chr)%>%distinct(kmermiddle, .keep_all=TRUE)%>%ungroup
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
df%>%arrange(-count)
gp1 <- ggplot(df%>%filter(annotation=="p_euchromatin")) + geom_rect(aes(xmin=kstart, xmax=kend, ymin=0, ymax=count, fill=strand_color))
gp1 <- gp1 + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.) + scale_alpha_identity() + theme_bw() 
gp1 <- gp1 + scale_x_continuous(breaks=c(-150000000, -100000000, -50000000,-5000000),limits=c(-148102972,-5000000)) + theme(strip.background = element_blank(), strip.text.y = element_blank(),legend.position = "none")
gp1 <- gp1 +scale_y_continuous(breaks = seq(-5000,10000, by=5000), limits=c(-5000,10000))
print(gp1)

gp2 <- ggplot(df%>%filter(annotation=="p_pericentromere")) + geom_rect(aes(xmin=kstart, xmax=kend, ymin=0, ymax=count, fill=strand_color))
gp2 <- gp2 + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.) + scale_alpha_identity() + theme_bw() 
gp2 <- gp2 + coord_cartesian(xlim=c(-5000000, 0)) + theme(strip.background = element_blank(), axis.text.y=element_blank(), strip.text.y = element_blank(),axis.ticks.y=element_blank(),legend.position = "none")
gp2 <- gp2 +scale_y_continuous(breaks = seq(-5000,10000, by=5000), limits=c(-5000,10000))

gp3 <- ggplot(df%>%filter(annotation=="q_pericentromere")) + geom_rect(aes(xmin=kstart, xmax=kend, ymin=0, ymax=count, fill=strand_color))
gp3 <- gp3 + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.) + scale_alpha_identity() + theme_bw() 
gp3 <- gp3 + coord_cartesian(xlim=c(0, 5000000)) + theme(strip.background = element_blank(), axis.text.y=element_blank(), strip.text.y = element_blank(),axis.ticks.y=element_blank(),legend.position = "none")
gp3 <- gp3 +scale_y_continuous(breaks = seq(-5000,10000, by=5000), limits=c(-5000,10000))

gp4 <- ggplot(df%>%filter(annotation=="q_pericentromere")) + geom_rect(aes(xmin=kstart, xmax=kend, ymin=0, ymax=count, fill=strand_color))
gp4 <- gp4 + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.) + scale_alpha_identity() + theme_bw() 
gp4 <- gp4 + scale_x_continuous(breaks=c(5000000, 50000000, 100000000, 150000000),limits=c(5000000,148102972)) + theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(), legend.position = "none")
gp4 <- gp4 +scale_y_continuous(breaks = seq(-5000,10000, by=5000), limits=c(-5000,10000))
print(gp4)

grid.arrange(gp1, gp2, ncol=2)

#grob1<-ggplotGrob(gp1)
#grob2<-ggplotGrob(gp2)
#grob3<-ggplotGrob(gp3)
#grob4<-ggplotGrob(gp4)

#grob1$widths<-grob4$widths*(1/5)
#grob2$widths<-grob4$widths*(1/2)
#grob3$widths<-grob4$widths*(1/2)
#g<-cbind(grob1, grob2, grob3, grob4)
g<-plot_grid(gp1, gp2, gp3, gp4, ncol=4, scale=1.01)
print(g)

(savefile = str_c(k, "mer_chr", hor, "HOR_hit_test.png",sep = ""))

ggsave(savefile, g, height=15, width=15)


k = 14
hor = 11
(gfffile = str_c(k, "mer_chr", hor, "HOR_hit_test.gff",sep = ""))
(kmerhit <- read_csv(gfffile, col_names = c("seq", "source", "feature", "kmerstart", "kmerend", "score", "strand", "attribute"))
  %>% filter(feature == "kmer_hit") %>%separate(seq, c("filename", "chr"), ":")
  %>% select(chr, kmerstart, kmerend, strand))
(kmerhit<-kmerhit%>%mutate("kmermiddle"=(kmerstart+kmerend)/2))

(kmer_rel<-kmerhit%>%left_join(cent)%>%filter(kmermiddle<CentStart)%>%mutate("kmer_rel"=kmermiddle-CentStart)
  %>%bind_rows(kmerhit%>%left_join(cent)%>%filter(kmermiddle>CentEnd)%>%mutate("kmer_rel"=kmermiddle-CentEnd)))
kmer_rel<-kmer_rel%>%group_by(chr)%>%distinct(kmermiddle, .keep_all=TRUE)%>%ungroup
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
df%>%arrange(-count)
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
gp1 <- gp1 + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.) + scale_alpha_identity() + theme_bw() 
gp1 <- gp1 + scale_x_continuous(breaks=c(-150000000, -100000000, -50000000,-5000000),limits=c(-148102972,-5000000)) + theme(strip.background = element_blank(), strip.text.y = element_blank(),legend.position = "none")
gp1 <- gp1 +scale_y_continuous(breaks = seq(-5000,15000, by=5000), limits=c(-5000,15000))
print(gp)

gp2 <- ggplot(df%>%filter(annotation=="p_pericentromere")) + geom_rect(aes(xmin=kstart, xmax=kend, ymin=0, ymax=count, fill=strand_color))
gp2 <-gp2 + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.) + scale_alpha_identity() + theme_bw() 
gp2 <-gp2 + coord_cartesian(xlim=c(-5000000, 0)) + theme(strip.background = element_blank(), axis.text.y=element_blank(), strip.text.y = element_blank(),axis.ticks.y=element_blank(),legend.position = "none")
gp2 <- gp2 +scale_y_continuous(breaks = seq(-5000,15000, by=5000), limits=c(-5000,15000))

gp3 <- ggplot(df%>%filter(annotation=="q_pericentromere")) + geom_rect(aes(xmin=kstart, xmax=kend, ymin=0, ymax=count, fill=strand_color))
gp3 <-gp3 + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.) + scale_alpha_identity() + theme_bw() 
gp3 <-gp3 + coord_cartesian(xlim=c(0, 5000000)) + theme(strip.background = element_blank(), axis.text.y=element_blank(), strip.text.y = element_blank(),axis.ticks.y=element_blank(),legend.position = "none")
gp3 <- gp3 +scale_y_continuous(breaks = seq(-5000,15000, by=5000), limits=c(-5000,15000))

gp4 <- ggplot(df%>%filter(annotation=="q_pericentromere")) + geom_rect(aes(xmin=kstart, xmax=kend, ymin=0, ymax=count, fill=strand_color))
gp4 <- gp4 + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.) + scale_alpha_identity() + theme_bw() 
gp4 <- gp4 + scale_x_continuous(breaks=c(5000000, 50000000, 100000000, 150000000),limits=c(5000000,148102972)) + theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(), legend.position = "none")
gp4 <- gp4 +scale_y_continuous(breaks = seq(-5000,15000, by=5000), limits=c(-5000,15000))


g<-plot_grid(gp1, gp2, gp3, gp4, ncol=4, scale=1.01)
print(g)

(savefile = str_c(k, "mer_chr", hor, "HOR_hit_test.png",sep = ""))

ggsave(savefile, g, height=15, width=15)



k = 18
hor = 1
(gfffile = str_c(k, "mer_chr", hor, "HOR_hit_test.gff",sep = ""))
(kmerhit <- read_csv(gfffile, col_names = c("seq", "source", "feature", "kmerstart", "kmerend", "score", "strand", "attribute"))
  %>% filter(feature == "kmer_hit") %>%separate(seq, c("filename", "chr"), ":")
  %>% select(chr, kmerstart, kmerend, strand))
(kmerhit<-kmerhit%>%mutate("kmermiddle"=(kmerstart+kmerend)/2))


(kmer_rel<-kmerhit%>%left_join(cent)%>%filter(kmermiddle<CentStart)%>%mutate("kmer_rel"=kmermiddle-CentStart)
  %>%bind_rows(kmerhit%>%left_join(cent)%>%filter(kmermiddle>CentEnd)%>%mutate("kmer_rel"=kmermiddle-CentEnd)))
kmer_rel<-kmer_rel%>%group_by(chr)%>%distinct(kmermiddle, .keep_all=TRUE)%>%ungroup
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
df%>%arrange(-count)
gp1 <- ggplot(df%>%filter(annotation=="p_euchromatin")) + geom_rect(aes(xmin=kstart, xmax=kend, ymin=0, ymax=count, fill=strand_color))
gp1 <- gp1 + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.) + scale_alpha_identity() + theme_bw() 
gp1 <- gp1 + scale_x_continuous(breaks=c(-150000000, -100000000, -50000000,-5000000),limits=c(-148102972,-5000000)) + theme(strip.background = element_blank(), strip.text.y = element_blank(),legend.position = "none")
gp1 <- gp1 +scale_y_continuous(breaks = seq(-5000,10000, by=5000), limits=c(-5000,10000))
print(gp1)

gp2 <- ggplot(df%>%filter(annotation=="p_pericentromere")) + geom_rect(aes(xmin=kstart, xmax=kend, ymin=0, ymax=count, fill=strand_color))
gp2 <- gp2 + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.) + scale_alpha_identity() + theme_bw() 
gp2 <- gp2 + coord_cartesian(xlim=c(-5000000, 0)) + theme(strip.background = element_blank(), axis.text.y=element_blank(), strip.text.y = element_blank(),axis.ticks.y=element_blank(),legend.position = "none")
gp2 <- gp2 +scale_y_continuous(breaks = seq(-5000,10000, by=5000), limits=c(-5000,10000))

gp3 <- ggplot(df%>%filter(annotation=="q_pericentromere")) + geom_rect(aes(xmin=kstart, xmax=kend, ymin=0, ymax=count, fill=strand_color))
gp3 <- gp3 + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.) + scale_alpha_identity() + theme_bw() 
gp3 <- gp3 + coord_cartesian(xlim=c(0, 5000000)) + theme(strip.background = element_blank(), axis.text.y=element_blank(), strip.text.y = element_blank(),axis.ticks.y=element_blank(),legend.position = "none")
gp3 <- gp3 +scale_y_continuous(breaks = seq(-5000,10000, by=5000), limits=c(-5000,10000))

gp4 <- ggplot(df%>%filter(annotation=="q_pericentromere")) + geom_rect(aes(xmin=kstart, xmax=kend, ymin=0, ymax=count, fill=strand_color))
gp4 <- gp4 + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.) + scale_alpha_identity() + theme_bw() 
gp4 <- gp4 + scale_x_continuous(breaks=c(5000000, 50000000, 100000000, 150000000),limits=c(5000000,148102972)) + theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(), legend.position = "none")
gp4 <- gp4 +scale_y_continuous(breaks = seq(-5000,10000, by=5000), limits=c(-5000,10000))
print(gp)

# grid.arrange(gp1, gp2, ncol=2)

#grob1<-ggplotGrob(gp1)
#grob2<-ggplotGrob(gp2)
#grob3<-ggplotGrob(gp3)
#grob4<-ggplotGrob(gp4)

#grob1$widths<-grob4$widths*(1/5)
#grob2$widths<-grob4$widths*(1/2)
#grob3$widths<-grob4$widths*(1/2)
#g<-cbind(grob1, grob2, grob3, grob4)
g<-plot_grid(gp1, gp2, gp3, gp4, ncol=4, scale=1.01)
print(g)

(savefile = str_c(k, "mer_chr", hor, "HOR_hit_test.png",sep = ""))

ggsave(savefile, g, height=15, width=15)


k = 14
hor = 1
(gfffile = str_c(k, "mer_chr", hor, "HOR_hit_test.gff",sep = ""))
(kmerhit <- read_csv(gfffile, col_names = c("seq", "source", "feature", "kmerstart", "kmerend", "score", "strand", "attribute"))
  %>% filter(feature == "kmer_hit") %>%separate(seq, c("filename", "chr"), ":")
  %>% select(chr, kmerstart, kmerend, strand))
(kmerhit<-kmerhit%>%mutate("kmermiddle"=(kmerstart+kmerend)/2))

(kmer_rel<-kmerhit%>%left_join(cent)%>%filter(kmermiddle<CentStart)%>%mutate("kmer_rel"=kmermiddle-CentStart)
  %>%bind_rows(kmerhit%>%left_join(cent)%>%filter(kmermiddle>CentEnd)%>%mutate("kmer_rel"=kmermiddle-CentEnd)))
kmer_rel<-kmer_rel%>%group_by(chr)%>%distinct(kmermiddle, .keep_all=TRUE)%>%ungroup
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
df%>%arrange(-count)
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
gp1 <- gp1 + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.) + scale_alpha_identity() + theme_bw() 
gp1 <- gp1 + scale_x_continuous(breaks=c(-150000000, -100000000, -50000000,-5000000),limits=c(-148102972,-5000000)) + theme(strip.background = element_blank(), strip.text.y = element_blank(),legend.position = "none")
gp1 <- gp1 +scale_y_continuous(breaks = seq(-10000,20000, by=10000), limits=c(-10000,20000))
print(gp)

gp2 <- ggplot(df%>%filter(annotation=="p_pericentromere")) + geom_rect(aes(xmin=kstart, xmax=kend, ymin=0, ymax=count, fill=strand_color))
gp2 <-gp2 + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.) + scale_alpha_identity() + theme_bw() 
gp2 <-gp2 + coord_cartesian(xlim=c(-5000000, 0)) + theme(strip.background = element_blank(), axis.text.y=element_blank(), strip.text.y = element_blank(),axis.ticks.y=element_blank(),legend.position = "none")
gp2 <- gp2 +scale_y_continuous(breaks = seq(-10000,20000, by=10000), limits=c(-10000,20000))

gp3 <- ggplot(df%>%filter(annotation=="q_pericentromere")) + geom_rect(aes(xmin=kstart, xmax=kend, ymin=0, ymax=count, fill=strand_color))
gp3 <-gp3 + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.) + scale_alpha_identity() + theme_bw() 
gp3 <-gp3 + coord_cartesian(xlim=c(0, 5000000)) + theme(strip.background = element_blank(), axis.text.y=element_blank(), strip.text.y = element_blank(),axis.ticks.y=element_blank(),legend.position = "none")
gp3 <- gp3 +scale_y_continuous(breaks = seq(-10000,20000, by=10000), limits=c(-10000,20000))

gp4 <- ggplot(df%>%filter(annotation=="q_pericentromere")) + geom_rect(aes(xmin=kstart, xmax=kend, ymin=0, ymax=count, fill=strand_color))
gp4 <- gp4 + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.) + scale_alpha_identity() + theme_bw() 
gp4 <- gp4 + scale_x_continuous(breaks=c(5000000, 50000000, 100000000, 150000000),limits=c(5000000,148102972)) + theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(), legend.position = "none")
gp4 <- gp4 +scale_y_continuous(breaks = seq(-10000,20000, by=10000), limits=c(-10000,20000))


g<-plot_grid(gp1, gp2, gp3, gp4, ncol=4, scale=1.01)
print(g)

(savefile = str_c(k, "mer_chr", hor, "HOR_hit_test.png",sep = ""))

ggsave(savefile, g, height=15, width=15)


k = 18
hor = 2
(gfffile = str_c(k, "mer_chr", hor, "HOR_hit_test.gff",sep = ""))
(kmerhit <- read_csv(gfffile, col_names = c("seq", "source", "feature", "kmerstart", "kmerend", "score", "strand", "attribute"))
  %>% filter(feature == "kmer_hit") %>%separate(seq, c("filename", "chr"), ":")
  %>% select(chr, kmerstart, kmerend, strand))
(kmerhit<-kmerhit%>%mutate("kmermiddle"=(kmerstart+kmerend)/2))


(kmer_rel<-kmerhit%>%left_join(cent)%>%filter(kmermiddle<CentStart)%>%mutate("kmer_rel"=kmermiddle-CentStart)
  %>%bind_rows(kmerhit%>%left_join(cent)%>%filter(kmermiddle>CentEnd)%>%mutate("kmer_rel"=kmermiddle-CentEnd)))
kmer_rel<-kmer_rel%>%group_by(chr)%>%distinct(kmermiddle, .keep_all=TRUE)%>%ungroup
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
df%>%arrange(-count)
gp1 <- ggplot(df%>%filter(annotation=="p_euchromatin")) + geom_rect(aes(xmin=kstart, xmax=kend, ymin=0, ymax=count, fill=strand_color))
gp1 <- gp1 + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.) + scale_alpha_identity() + theme_bw() 
gp1 <- gp1 + scale_x_continuous(breaks=c(-150000000, -100000000, -50000000,-5000000),limits=c(-148102972,-5000000)) + theme(strip.background = element_blank(), strip.text.y = element_blank(),legend.position = "none")
gp1 <- gp1 +scale_y_continuous(breaks = seq(-10000,20000, by=10000), limits=c(-10000,20000))
print(gp)

gp2 <- ggplot(df%>%filter(annotation=="p_pericentromere")) + geom_rect(aes(xmin=kstart, xmax=kend, ymin=0, ymax=count, fill=strand_color))
gp2 <- gp2 + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.) + scale_alpha_identity() + theme_bw() 
gp2 <- gp2 + coord_cartesian(xlim=c(-5000000, 0)) + theme(strip.background = element_blank(), axis.text.y=element_blank(), strip.text.y = element_blank(),axis.ticks.y=element_blank(),legend.position = "none")
gp2 <- gp2 +scale_y_continuous(breaks = seq(-10000,20000, by=10000), limits=c(-10000,20000))

gp3 <- ggplot(df%>%filter(annotation=="q_pericentromere")) + geom_rect(aes(xmin=kstart, xmax=kend, ymin=0, ymax=count, fill=strand_color))
gp3 <- gp3 + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.) + scale_alpha_identity() + theme_bw() 
gp3 <- gp3 + coord_cartesian(xlim=c(0, 5000000)) + theme(strip.background = element_blank(), axis.text.y=element_blank(), strip.text.y = element_blank(),axis.ticks.y=element_blank(),legend.position = "none")
gp3 <- gp3 +scale_y_continuous(breaks = seq(-10000,20000, by=10000), limits=c(-10000,20000))

gp4 <- ggplot(df%>%filter(annotation=="q_pericentromere")) + geom_rect(aes(xmin=kstart, xmax=kend, ymin=0, ymax=count, fill=strand_color))
gp4 <- gp4 + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.) + scale_alpha_identity() + theme_bw() 
gp4 <- gp4 + scale_x_continuous(breaks=c(5000000, 50000000, 100000000, 150000000),limits=c(5000000,148102972)) + theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(), legend.position = "none")
gp4 <- gp4 +scale_y_continuous(breaks = seq(-10000,20000, by=10000), limits=c(-10000,20000))
print(gp)

# grid.arrange(gp1, gp2, ncol=2)

#grob1<-ggplotGrob(gp1)
#grob2<-ggplotGrob(gp2)
#grob3<-ggplotGrob(gp3)
#grob4<-ggplotGrob(gp4)

#grob1$widths<-grob4$widths*(1/5)
#grob2$widths<-grob4$widths*(1/2)
#grob3$widths<-grob4$widths*(1/2)
#g<-cbind(grob1, grob2, grob3, grob4)
g<-plot_grid(gp1, gp2, gp3, gp4, ncol=4, scale=1.01)
print(g)

(savefile = str_c(k, "mer_chr", hor, "HOR_hit_test.png",sep = ""))

ggsave(savefile, g, height=15, width=15)


k = 14
hor = 2
(gfffile = str_c(k, "mer_chr", hor, "HOR_hit_test.gff",sep = ""))
(kmerhit <- read_csv(gfffile, col_names = c("seq", "source", "feature", "kmerstart", "kmerend", "score", "strand", "attribute"))
  %>% filter(feature == "kmer_hit") %>%separate(seq, c("filename", "chr"), ":")
  %>% select(chr, kmerstart, kmerend, strand))
(kmerhit<-kmerhit%>%mutate("kmermiddle"=(kmerstart+kmerend)/2))

(kmer_rel<-kmerhit%>%left_join(cent)%>%filter(kmermiddle<CentStart)%>%mutate("kmer_rel"=kmermiddle-CentStart)
  %>%bind_rows(kmerhit%>%left_join(cent)%>%filter(kmermiddle>CentEnd)%>%mutate("kmer_rel"=kmermiddle-CentEnd)))
kmer_rel<-kmer_rel%>%group_by(chr)%>%distinct(kmermiddle, .keep_all=TRUE)%>%ungroup
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
df%>%arrange(-count)
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
gp1 <- gp1 + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.) + scale_alpha_identity() + theme_bw() 
gp1 <- gp1 + scale_x_continuous(breaks=c(-150000000, -100000000, -50000000,-5000000),limits=c(-148102972,-5000000)) + theme(strip.background = element_blank(), strip.text.y = element_blank(),legend.position = "none")
gp1 <- gp1 +scale_y_continuous(breaks = seq(-10000,20000, by=10000), limits=c(-10000,20000))
print(gp)

gp2 <- ggplot(df%>%filter(annotation=="p_pericentromere")) + geom_rect(aes(xmin=kstart, xmax=kend, ymin=0, ymax=count, fill=strand_color))
gp2 <-gp2 + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.) + scale_alpha_identity() + theme_bw() 
gp2 <-gp2 + coord_cartesian(xlim=c(-5000000, 0)) + theme(strip.background = element_blank(), axis.text.y=element_blank(), strip.text.y = element_blank(),axis.ticks.y=element_blank(),legend.position = "none")
gp2 <- gp2 +scale_y_continuous(breaks = seq(-10000,20000, by=10000), limits=c(-10000,20000))

gp3 <- ggplot(df%>%filter(annotation=="q_pericentromere")) + geom_rect(aes(xmin=kstart, xmax=kend, ymin=0, ymax=count, fill=strand_color))
gp3 <-gp3 + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.) + scale_alpha_identity() + theme_bw() 
gp3 <-gp3 + coord_cartesian(xlim=c(0, 5000000)) + theme(strip.background = element_blank(), axis.text.y=element_blank(), strip.text.y = element_blank(),axis.ticks.y=element_blank(),legend.position = "none")
gp3 <- gp3 +scale_y_continuous(breaks = seq(-10000,20000, by=10000), limits=c(-10000,20000))

gp4 <- ggplot(df%>%filter(annotation=="q_pericentromere")) + geom_rect(aes(xmin=kstart, xmax=kend, ymin=0, ymax=count, fill=strand_color))
gp4 <- gp4 + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.) + scale_alpha_identity() + theme_bw() 
gp4 <- gp4 + scale_x_continuous(breaks=c(5000000, 50000000, 100000000, 150000000),limits=c(5000000,148102972)) + theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(), legend.position = "none")
gp4 <- gp4 +scale_y_continuous(breaks = seq(-10000,20000, by=10000), limits=c(-10000,20000))


g<-plot_grid(gp1, gp2, gp3, gp4, ncol=4, scale=1.01)
print(g)

(savefile = str_c(k, "mer_chr", hor, "HOR_hit_test.png",sep = ""))

ggsave(savefile, g, height=15, width=15)






k = 19
hor = 1
(kmerfile = str_c("chr", hor, "_", k, "merlist.csv",sep = ""))
(kmerlist <- read_csv(kmerfile, col_names = c("index", "seq")))
for (hor in c(2, 11)){
  (kmerfile = str_c("chr", hor, "_", k, "merlist.csv",sep = ""))
  (l<-read_csv(kmerfile, col_names = c("index", "seq")))
  (kmerlist <- kmerlist%>%full_join(l, by="seq"))
}
kmerlist%>%filter(is.na(index.x)==FALSE & is.na(index)==FALSE|is.na(index.y)==FALSE & is.na(index)==FALSE|is.na(index.x)==FALSE & is.na(index.y)==FALSE)%>%View

# kmerlist%>%filter(is.na(index.x) & is.na(index)|is.na(index.y) & is.na(index)|is.na(index.x) & is.na(index.y))%>%View
kmerlist1<-kmerlist%>%filter(is.na(index.y) & is.na(index))%>%rename("attribute" = index.x)%>%select(seq,attribute)
kmerlist2<-kmerlist%>%filter(is.na(index.x) & is.na(index))%>%rename("attribute" = index.y)%>%select(seq,attribute)
kmerlist11<-kmerlist%>%filter(is.na(index.x) & is.na(index.y))%>%rename("attribute" = index)%>%select(seq,attribute)

594/600
882/899
801/834

(peri<-read_csv("peri_and_livingcent.bed", col_names = c("chr", "start", "end")))
(cent<-read_csv("cent_annotation.csv")%>%filter(State=="live")%>%select(chr, CentStart, CentEnd))

hor = 1
(gfffile = str_c(k, "mer_chr", hor, "HOR_hit_test.gff",sep = ""))
(kmerhit <- read_csv(gfffile, col_names = c("seq", "source", "feature", "kmerstart", "kmerend", "score", "strand", "attribute"))
  %>% semi_join(., kmerlist1, by="attribute")%>%group_by(seq)%>%distinct(kmerstart, .keep_all = TRUE)%>%ungroup()
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

(kmer_rel<-kmerhit%>%left_join(cent)%>%filter(kmermiddle<CentStart)%>%mutate("kmer_rel"=kmermiddle-CentStart)
  %>%bind_rows(kmerhit%>%left_join(cent)%>%filter(kmermiddle>CentEnd)%>%mutate("kmer_rel"=kmermiddle-CentEnd)))
kmer_rel<-kmer_rel%>%group_by(chr)%>%distinct(kmermiddle, .keep_all=TRUE)%>%ungroup
(cnt_rel<-kmer_rel%>%filter(kmer_rel>0)%>%filter(strand=="+")%>%group_by(chr)%>%count(cut_width(kmer_rel, 100000, boundary=0))%>%ungroup%>%mutate(strand=1, strand_color="+")%>%
    bind_rows(kmer_rel%>%filter(kmer_rel>0)%>%filter(strand=="-")%>%group_by(chr)%>%count(cut_width(kmer_rel, 100000, boundary=0))%>%ungroup%>%mutate(strand=-1, strand_color="-"))%>%
    bind_rows(kmer_rel%>%filter(kmer_rel<0)%>%filter(strand=="+")%>%group_by(chr)%>%count(cut_width(kmer_rel, 100000, boundary=0))%>%ungroup%>%mutate(strand=1, strand_color="+"))%>%
    bind_rows(kmer_rel%>%filter(kmer_rel<0)%>%filter(strand=="-")%>%group_by(chr)%>%count(cut_width(kmer_rel, 100000, boundary=0))%>%ungroup%>%mutate(strand=-1, strand_color="-"))%>%
    mutate("position"=str_sub(`cut_width(kmer_rel, 1e+05, boundary = 0)`, 2, -2))%>%separate(position,c("kstart", "kend"),",")%>%
    mutate("count"=n*strand))


cnt_rel$kstart<-as.integer(cnt_rel$kstart)
cnt_rel$kend<-as.integer(cnt_rel$kend)
cnt_rel<-cnt_rel%>%mutate("x_rel"=(kstart+kend)/2)
cnt_rel1<-cnt_rel%>%mutate("SF_color"="SF1")

hor = 2
(gfffile = str_c(k, "mer_chr", hor, "HOR_hit_test.gff",sep = ""))
(kmerhit <- read_csv(gfffile, col_names = c("seq", "source", "feature", "kmerstart", "kmerend", "score", "strand", "attribute"))
  %>% semi_join(., kmerlist2, by="attribute")%>%group_by(seq)%>%distinct(kmerstart, .keep_all = TRUE)%>%ungroup()
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

(kmer_rel<-kmerhit%>%left_join(cent)%>%filter(kmermiddle<CentStart)%>%mutate("kmer_rel"=kmermiddle-CentStart)
  %>%bind_rows(kmerhit%>%left_join(cent)%>%filter(kmermiddle>CentEnd)%>%mutate("kmer_rel"=kmermiddle-CentEnd)))
kmer_rel<-kmer_rel%>%group_by(chr)%>%distinct(kmermiddle, .keep_all=TRUE)%>%ungroup
(cnt_rel<-kmer_rel%>%filter(kmer_rel>0)%>%filter(strand=="+")%>%group_by(chr)%>%count(cut_width(kmer_rel, 100000, boundary=0))%>%ungroup%>%mutate(strand=1, strand_color="+")%>%
    bind_rows(kmer_rel%>%filter(kmer_rel>0)%>%filter(strand=="-")%>%group_by(chr)%>%count(cut_width(kmer_rel, 100000, boundary=0))%>%ungroup%>%mutate(strand=-1, strand_color="-"))%>%
    bind_rows(kmer_rel%>%filter(kmer_rel<0)%>%filter(strand=="+")%>%group_by(chr)%>%count(cut_width(kmer_rel, 100000, boundary=0))%>%ungroup%>%mutate(strand=1, strand_color="+"))%>%
    bind_rows(kmer_rel%>%filter(kmer_rel<0)%>%filter(strand=="-")%>%group_by(chr)%>%count(cut_width(kmer_rel, 100000, boundary=0))%>%ungroup%>%mutate(strand=-1, strand_color="-"))%>%
    mutate("position"=str_sub(`cut_width(kmer_rel, 1e+05, boundary = 0)`, 2, -2))%>%separate(position,c("kstart", "kend"),",")%>%
    mutate("count"=n*strand))


cnt_rel$kstart<-as.integer(cnt_rel$kstart)
cnt_rel$kend<-as.integer(cnt_rel$kend)
cnt_rel<-cnt_rel%>%mutate("x_rel"=(kstart+kend)/2)
cnt_rel2<-cnt_rel%>%mutate("SF_color"="SF2")


hor = 11
(gfffile = str_c(k, "mer_chr", hor, "HOR_hit_test.gff",sep = ""))
(kmerhit <- read_csv(gfffile, col_names = c("seq", "source", "feature", "kmerstart", "kmerend", "score", "strand", "attribute"))
  %>% semi_join(., kmerlist11, by="attribute")%>%group_by(seq)%>%distinct(kmerstart, .keep_all = TRUE)%>%ungroup()
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

(kmer_rel<-kmerhit%>%left_join(cent)%>%filter(kmermiddle<CentStart)%>%mutate("kmer_rel"=kmermiddle-CentStart)
  %>%bind_rows(kmerhit%>%left_join(cent)%>%filter(kmermiddle>CentEnd)%>%mutate("kmer_rel"=kmermiddle-CentEnd)))
kmer_rel<-kmer_rel%>%group_by(chr)%>%distinct(kmermiddle, .keep_all=TRUE)%>%ungroup
(cnt_rel<-kmer_rel%>%filter(kmer_rel>0)%>%filter(strand=="+")%>%group_by(chr)%>%count(cut_width(kmer_rel, 100000, boundary=0))%>%ungroup%>%mutate(strand=1, strand_color="+")%>%
    bind_rows(kmer_rel%>%filter(kmer_rel>0)%>%filter(strand=="-")%>%group_by(chr)%>%count(cut_width(kmer_rel, 100000, boundary=0))%>%ungroup%>%mutate(strand=-1, strand_color="-"))%>%
    bind_rows(kmer_rel%>%filter(kmer_rel<0)%>%filter(strand=="+")%>%group_by(chr)%>%count(cut_width(kmer_rel, 100000, boundary=0))%>%ungroup%>%mutate(strand=1, strand_color="+"))%>%
    bind_rows(kmer_rel%>%filter(kmer_rel<0)%>%filter(strand=="-")%>%group_by(chr)%>%count(cut_width(kmer_rel, 100000, boundary=0))%>%ungroup%>%mutate(strand=-1, strand_color="-"))%>%
    mutate("position"=str_sub(`cut_width(kmer_rel, 1e+05, boundary = 0)`, 2, -2))%>%separate(position,c("kstart", "kend"),",")%>%
    mutate("count"=n*strand))


cnt_rel$kstart<-as.integer(cnt_rel$kstart)
cnt_rel$kend<-as.integer(cnt_rel$kend)
cnt_rel<-cnt_rel%>%mutate("x_rel"=(kstart+kend)/2)
cnt_rel11<-cnt_rel%>%mutate("SF_color"="SF3")


(df<-bind_rows(cnt_rel1, cnt_rel2, cnt_rel11))
df$chr_f = factor(df$chr, levels=chromosomes)

df%>%filter(chr=="chr1"|chr=="chr3"|chr=="chr5"|chr=="chr6"|chr=="chr7"|chr=="chr10"|chr=="chr12"|chr=="chr16"|chr=="chr19")%>%mutate(sum(count))
df%>%filter(chr=="chr1"|chr=="chr3"|chr=="chr5"|chr=="chr6"|chr=="chr7"|chr=="chr10"|chr=="chr12"|chr=="chr16"|chr=="chr19")%>%filter(SF_color=="SF1")%>%mutate(sum(count))
78971/198710

df%>%filter(chr=="chr2"|chr=="chr4"|chr=="chr8"|chr=="chr9"|chr=="chr13"|chr=="chr14"|chr=="chr15"|chr=="chr18"|chr=="chr20"|chr=="chr21"|chr=="chr22")%>%mutate(sum(count))
df%>%filter(chr=="chr2"|chr=="chr4"|chr=="chr8"|chr=="chr9"|chr=="chr13"|chr=="chr14"|chr=="chr15"|chr=="chr18"|chr=="chr20"|chr=="chr21"|chr=="chr22")%>%filter(SF_color=="SF2")%>%mutate(sum(count))
115261/163252

df%>%filter(chr=="chr11"|chr=="chr17"|chr=="chrX")%>%mutate(sum(count))
df%>%filter(chr=="chr11"|chr=="chr17"|chr=="chrX")%>%filter(SF_color=="SF3")%>%mutate(sum(count))
29112/51661



df%>%filter(SF_color=="SF1")%>%mutate(sum(count))
df%>%filter(SF_color=="SF1")%>%filter(chr=="chr1"|chr=="chr3"|chr=="chr5"|chr=="chr6"|chr=="chr7"|chr=="chr10"|chr=="chr12"|chr=="chr16"|chr=="chr19")%>%mutate(sum(count))

df%>%filter(SF_color=="SF2")%>%mutate(sum(count))
df%>%filter(SF_color=="SF2")%>%filter(chr=="chr2"|chr=="chr4"|chr=="chr8"|chr=="chr9"|chr=="chr13"|chr=="chr14"|chr=="chr15"|chr=="chr18"|chr=="chr20"|chr=="chr21"|chr=="chr22")%>%mutate(sum(count))

df%>%filter(SF_color=="SF3")%>%mutate(sum(count))
df%>%filter(SF_color=="SF3")%>%filter(chr=="chr11"|chr=="chr17"|chr=="chrX")%>%mutate(sum(count))

78981/89472
115261/199463
29112/121777
df1<-df%>%filter(x_rel<0)%>%mutate("PQ"="p-arm")%>%bind_rows(df%>%filter(x_rel>=0)%>%mutate("PQ"="q-arm"))

chromosomes<-c("chr1","chr3","chr5","chr6","chr7","chr10","chr12","chr16","chr19","chr2","chr4", "chr8", "chr9","chr13","chr14","chr15","chr18","chr20","chr21","chr22", "chr11", "chr17", "chrX", "chrY")
df1$chr_f = factor(df1$chr, levels=chromosomes)


gp <-  ggplot(df1%>%filter(x_rel>-3000000 & x_rel<3000000)) + geom_col(aes(x=x_rel, y=count, fill=SF_color), position="stack")
gp<- gp + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~PQ, scales = "free_x")+xlab("relative distance") # +coord_cartesian(xlim=c(-3000000, 3000000))
# gp<- gp + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.)+coord_cartesian(xlim=c(-3000000, 0))+xlab("relative distance")
gp <- gp + scale_alpha_identity() + theme_bw() +scale_fill_hue(labels=c(SF1="chr1", SF2="chr2", SF3="chr11"))+  theme(legend.title=element_blank())
print(gp)
ggsave("SF_kmerhit.png", gp, width=8, height=9)



gp <-  ggplot(df1%>%filter(x_rel>-3000000 & x_rel<3000000)%>%filter(chr=="chr1"|chr=="chr3"|chr=="chr5"|chr=="chr6"|chr=="chr7"|chr=="chr10"|chr=="chr12"|chr=="chr16"|chr=="chr19")) + geom_col(aes(x=x_rel, y=count, fill=SF_color), position="stack")
gp<- gp + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~PQ, scales = "free_x")+xlab("relative distance") # +coord_cartesian(xlim=c(-3000000, 3000000))
# gp<- gp + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.)+coord_cartesian(xlim=c(-3000000, 0))+xlab("relative distance")
gp <- gp + scale_alpha_identity() + theme_bw() +scale_fill_hue(labels=c(SF1="chr1", SF2="chr2", SF3="chr11"))+  theme(legend.title=element_blank())
print(gp)
ggsave("SF_kmerhit_SF1.png", gp, width=8, height=9)


gp <-  ggplot(df1%>%filter(x_rel>-3000000 & x_rel<3000000)%>%filter(chr=="chr2"|chr=="chr4"|chr=="chr8"|chr=="chr9"|chr=="chr13"|chr=="chr14"|chr=="chr15"|chr=="chr18"|chr=="chr20"|chr=="chr21"|chr=="chr22")) + geom_col(aes(x=x_rel, y=count, fill=SF_color), position="stack")
gp<- gp + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~PQ, scales = "free_x")+xlab("relative distance") # +coord_cartesian(xlim=c(-3000000, 3000000))
# gp<- gp + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.)+coord_cartesian(xlim=c(-3000000, 0))+xlab("relative distance")
gp <- gp + scale_alpha_identity() + theme_bw() +scale_fill_hue(labels=c(SF1="chr1", SF2="chr2", SF3="chr11"))+  theme(legend.title=element_blank())
print(gp)
ggsave("SF_kmerhit_SF2.png", gp, width=8, height=9)


gp <-  ggplot(df1%>%filter(x_rel>-3000000 & x_rel<3000000)%>%filter(chr=="chr11"|chr=="chr17"|chr=="chrX")) + geom_col(aes(x=x_rel, y=count, fill=SF_color), position="stack")
gp<- gp + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~PQ, scales = "free_x")+xlab("relative distance") # +coord_cartesian(xlim=c(-3000000, 3000000))
# gp<- gp + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.)+coord_cartesian(xlim=c(-3000000, 0))+xlab("relative distance")
gp <- gp + scale_alpha_identity() + theme_bw() +scale_fill_hue(labels=c(SF1="chr1", SF2="chr2", SF3="chr11"))+  theme(legend.title=element_blank())
print(gp)
ggsave("SF_kmerhit_SF3.png", gp, width=8, height=9)














gp <-  ggplot(df1%>%filter(SF_color=="SF1")%>%filter(x_rel>-3000000 & x_rel<3000000)) + geom_col(aes(x=x_rel, y=count, fill=SF_color), position="stack")+ scale_fill_manual(values=c("#F8766D"),labels=c(SF1="chr1"))
gp<- gp + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~PQ, scales = "free_x")+xlab("relative distance") #+coord_cartesian(xlim=c(-3000000, 3000000))
gp <- gp + scale_alpha_identity() + theme_bw() + theme(legend.title=element_blank())
print(gp)
ggsave("SF1_kmerhit.png", gp, width=7, height=10)

gp <-  ggplot(df1%>%filter(SF_color=="SF2")%>%filter(x_rel>-3000000 & x_rel<3000000)) + geom_col(aes(x=x_rel, y=count, fill=SF_color), position="stack") + scale_fill_manual(values=c("#7CAE00"), labels=c(SF2="chr2"))
gp<- gp + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~PQ, scales = "free_x")+xlab("relative distance") #+coord_cartesian(xlim=c(-3000000, 3000000))
gp <- gp + scale_alpha_identity() + theme_bw() + theme(legend.title=element_blank())
print(gp)
ggsave("SF2_kmerhit.png", gp,  width=7, height=10)

gp <-  ggplot(df1%>%filter(SF_color=="SF3")%>%filter(x_rel>-3000000 & x_rel<3000000)) + geom_col(aes(x=x_rel, y=count, fill=SF_color), position="stack") + scale_fill_manual(values=c("#00BFC4"),labels=c(SF3="chr11"))
gp<- gp + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~PQ, scales = "free_x")+xlab("relative distance") #+coord_cartesian(xlim=c(-3000000, 3000000))
gp <- gp + scale_alpha_identity() + theme_bw() + theme(legend.title=element_blank())
print(gp)
ggsave("SF3_kmerhit.png", gp,  width=7, height=10)



# gp <-  ggplot(df) + geom_rect(aes(xmin=kstart, xmax=kend, ymin=0, ymax=count, fill=SF_color))
gp <-  ggplot(df%>%filter(chr=="chr1"|chr=="chr3"|chr=="chr5"|chr=="chr6"|chr=="chr7"|chr=="chr10"|chr=="chr12"|chr=="chr16"|chr=="chr19")) + geom_col(aes(x=x_rel, y=count, fill=SF_color), position="stack")
gp<- gp + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.)+coord_cartesian(xlim=c(-5000000, 5000000))+xlab("relative distance")
gp <- gp + scale_alpha_identity() + theme_bw() 
print(gp)
ggsave("SF_kmerhit_SF1.png", gp, width=5, height=5)

gp <-  ggplot(df%>%filter(chr=="chr2"|chr=="chr4"|chr=="chr8"|chr=="chr9"|chr=="chr13"|chr=="chr14"|chr=="chr15"|chr=="chr18"|chr=="chr20"|chr=="chr21"|chr=="chr22")) + geom_col(aes(x=x_rel, y=count, fill=SF_color), position="stack")
gp<- gp + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.)+coord_cartesian(xlim=c(-5000000, 5000000))+xlab("relative distance")
gp <- gp + scale_alpha_identity() + theme_bw() 
print(gp)
ggsave("SF_kmerhit_SF2.png", gp, width=5, height=5)

gp <-  ggplot(df%>%filter(chr=="chr11"|chr=="chr17"|chr=="chrX")) + geom_col(aes(x=x_rel, y=count, fill=SF_color), position="stack")
gp<- gp + geom_vline(aes(xintercept=0, alpha=0.5), colour="#999999") + facet_grid(chr_f~.)+coord_cartesian(xlim=c(-5000000, 5000000))+xlab("relative distance")
gp <- gp + scale_alpha_identity() + theme_bw() 
print(gp)
ggsave("SF_kmerhit_SF3.png", gp, width=5, height=4)

