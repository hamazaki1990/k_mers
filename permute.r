library(tidyverse)
setwd("~/github/k_mers/")

hor = 1
k = 10
listfile = str_c(k, "mer_chr", hor, "HOR_randomcount.csv",sep = "")
origin = str_c("chr", hor,sep = "")
rowdata <-read_csv(listfile, col_names = c("EO", "cnt"))%>%mutate("chr" = origin, "kmer"=k)
rowdata

for (i in 11:21){
  k = i
  listfile = str_c(k, "mer_chr", hor, "HOR_randomcount.csv",sep = "")
  origin = str_c("chr", hor,sep = "")
  rowdata<-rowdata%>%bind_rows(read_csv(listfile, col_names = c("EO", "cnt"))%>%mutate("chr" = origin, "kmer"=k))
}
observed<-rowdata%>%filter(EO=="kmer")
df<-rowdata%>%group_by(chr)%>%group_by(kmer)%>%filter(EO!="kmer")%>%summarise("cnt"=mean(cnt))%>%mutate("chr"=origin)

for (j in c(2, 11)){
  hor = j
  k = 10
  listfile = str_c(k, "mer_chr", hor, "HOR_randomcount.csv",sep = "")
  origin = str_c("chr", hor,sep = "")
  rowdata <-read_csv(listfile, col_names = c("EO", "cnt"))%>%mutate("chr" = origin, "kmer"=k)
  rowdata

  for (i in 11:21){
    k = i
    listfile = str_c(k, "mer_chr", hor, "HOR_randomcount.csv",sep = "")
    origin = str_c("chr", hor,sep = "")
    rowdata<-rowdata%>%bind_rows(read_csv(listfile, col_names = c("EO", "cnt"))%>%mutate("chr" = origin, "kmer"=k))
  }
  observed<-observed%>%bind_rows(rowdata%>%filter(EO=="kmer"))
  df<-df%>%bind_rows(rowdata%>%group_by(chr)%>%group_by(kmer)%>%filter(EO!="kmer")%>%summarise("cnt"=mean(cnt))%>%mutate("chr"=origin))
}

#rowdata%>%View
df<-df%>%mutate("EO"="random")
df<-df%>%bind_rows(observed)
df%>%filter(EO=="kmer")
df$chr_f = factor(df$chr, levels=c("chr1", "chr2", "chr11"))
gp<-ggplot(df)+ geom_line(aes(x=kmer, y=cnt, color=EO))+facet_grid(chr_f~., scales="free")+theme_classic()
gp <- gp+scale_color_hue(labels=c(kmer="Observed", random="Expected(random)"))+  theme(legend.title=element_blank())
print(gp)
ggsave("kmerhit_EO.png", gp, width=10, height=10)
# df%>%View

gp<-ggplot(df%>%filter(chr=="chr1"))+ geom_line(aes(x=kmer, y=cnt, color=EO))+theme_classic()
gp <- gp+scale_color_hue(labels=c(kmer="Observed", random="Expected(random)"))+  theme(legend.title=element_blank())+xlab("k")+ylab("hit")
print(gp)
ggsave("kmerhit_EO_chr1.png", gp, width=4, height=3)



hor = 1
k = 10
listfile = str_c(k, "mer_chr", hor, "HOR_random_in_peri.csv",sep = "")
origin = str_c("chr", hor,sep = "")
rowdata_p <-read_csv(listfile, col_names = c("EO", "cnt"))%>%mutate("chr" = origin, "kmer"=k)
rowdata_p

for (i in 11:21){
  k = i
  listfile = str_c(k, "mer_chr", hor, "HOR_random_in_peri.csv",sep = "")
  origin = str_c("chr", hor,sep = "")
  rowdata_p<-rowdata_p%>%bind_rows(read_csv(listfile, col_names = c("EO", "cnt"))%>%mutate("chr" = origin, "kmer"=k))
}
observed<-rowdata_p%>%filter(EO=="kmer")
df_p<-rowdata_p%>%group_by(chr)%>%group_by(kmer)%>%filter(EO!="kmer")%>%summarise("cnt"=mean(cnt))%>%mutate("chr"=origin)

for (j in c(2, 11)){
  hor = j
  k = 10
  listfile = str_c(k, "mer_chr", hor, "HOR_random_in_peri.csv",sep = "")
  origin = str_c("chr", hor,sep = "")
  rowdata_p <-read_csv(listfile, col_names = c("EO", "cnt"))%>%mutate("chr" = origin, "kmer"=k)
  rowdata_p
  
  for (i in 11:21){
    k = i
    listfile = str_c(k, "mer_chr", hor, "HOR_random_in_peri.csv",sep = "")
    origin = str_c("chr", hor,sep = "")
    rowdata_p<-rowdata_p%>%bind_rows(read_csv(listfile, col_names = c("EO", "cnt"))%>%mutate("chr" = origin, "kmer"=k))
  }
  observed<-observed%>%bind_rows(rowdata_p%>%filter(EO=="kmer"))
  df_p<-df_p%>%bind_rows(rowdata_p%>%group_by(chr)%>%group_by(kmer)%>%filter(EO!="kmer")%>%summarise("cnt"=mean(cnt))%>%mutate("chr"=origin))
}

#rowdata_p%>%View
df_p<-df_p%>%mutate("EO"="random")
df_p<-df_p%>%bind_rows(observed)
df_p%>%filter(EO=="kmer")
df_p$chr_f = factor(df_p$chr, levels=c("chr1", "chr2", "chr11"))
(df_all<-df%>%filter(chr=="chr1")%>%filter(EO=="random")%>%left_join(df_p%>%filter(chr=="chr1")%>%filter(EO=="random")%>%select(cnt, kmer), by="kmer")%>%mutate("kmerinreg"=cnt.x-cnt.y)
  %>%bind_rows(df%>%filter(chr=="chr1")%>%filter(EO=="kmer")%>%left_join(df_p%>%filter(chr=="chr1")%>%filter(EO=="kmer")%>%select(cnt, kmer), by="kmer")%>%mutate("kmerinreg"=cnt.x-cnt.y)))
gp<-ggplot(df_p)+ geom_line(aes(x=kmer, y=cnt, color=EO))+facet_grid(chr_f~., scales="free")+theme_classic()
gp <- gp+scale_color_hue(labels=c(kmer="Observed", random="Expected(random)"))+  theme(legend.title=element_blank())
print(gp)
ggsave("kmerhit_EO_in_peri.png", gp, width=10, height=10)
# df_p%>%View
# df_all%>%gather(f_key, hit, cnt.y, kmerinreg)%>%View

gp<-ggplot(df_all%>%gather(f_key, hit, cnt.y, kmerinreg))+ geom_col(aes(x=kmer, y=hit, fill=EO), position ="dodge")+theme_classic()+facet_grid(.~f_key, labeller=as_labeller(c("cnt.y"="pericentromere", "kmerinreg"="out of pericentromere")))
gp <- gp+scale_fill_hue(labels=c(kmer="Observed", random="Expected"))+  theme(legend.title=element_blank()) + xlab("k") +ylab("hit")
print(gp)

ggsave("kmerhit_EO_out_peri_chr1_bar.png", gp, width=8, height=5)

gp<-ggplot(df_p%>%filter(chr=="chr1"))+ geom_col(aes(x=kmer, y=cnt, color=EO))+theme_classic()
gp <- gp+scale_color_hue(labels=c(kmer="Observed", random="Expected(random)"))+  theme(legend.title=element_blank())+xlab("k")+ylab("hit")
print(gp)
ggsave("kmerhit_EO_in_peri_chr1.png", gp, width=4, height=3)


(df_all%>%filter(kmer==19))

(9248550-6914914)/sqrt(6914914/10)

(2777031-1818158)/sqrt(1818158/10)

(1005315-484046)/sqrt(484046/10)

