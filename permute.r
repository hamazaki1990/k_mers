library(tidyverse)
setwd("~/github/k_mers/")

hor = 1
k = 9
listfile = str_c(k, "mer_chr", hor, "HOR_randomcount.csv",sep = "")
origin = str_c("chr", hor,sep = "")
rowdata <-read_csv(listfile, col_names = c("EO", "cnt"))%>%mutate("chr" = origin, "kmer"=k)
rowdata

for (i in 10:21){
  k = i
  listfile = str_c(k, "mer_chr", hor, "HOR_randomcount.csv",sep = "")
  origin = str_c("chr", hor,sep = "")
  rowdata<-rowdata%>%bind_rows(read_csv(listfile, col_names = c("EO", "cnt"))%>%mutate("chr" = origin, "kmer"=k))
}
observed<-rowdata%>%filter(EO=="kmer")
df<-rowdata%>%group_by(chr)%>%group_by(kmer)%>%filter(EO!="kmer")%>%summarise("cnt"=mean(cnt))%>%mutate("chr"=origin)

for (j in c(2, 11)){
  hor = j
  k = 9
  listfile = str_c(k, "mer_chr", hor, "HOR_randomcount.csv",sep = "")
  origin = str_c("chr", hor,sep = "")
  rowdata <-read_csv(listfile, col_names = c("EO", "cnt"))%>%mutate("chr" = origin, "kmer"=k)
  rowdata

  for (i in 10:21){
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
gp<-ggplot(df)+ geom_line(aes(x=kmer, y=cnt, color=EO))+facet_grid(chr~., scales="free")
print(gp)
ggsave("kmerhit_EO.png", gp, width=10, height=10)
# df%>%View


hor = 1
k = 15
listfile = str_c(k, "mer_chr", hor, "HOR_randomcount.csv",sep = "")
origin = str_c("chr", hor,sep = "")
rowdata <-read_csv(listfile, col_names = c("EO", "cnt"))%>%mutate("chr" = origin, "kmer"=k)
rowdata

for (i in 16:21){
  k = i
  listfile = str_c(k, "mer_chr", hor, "HOR_randomcount.csv",sep = "")
  origin = str_c("chr", hor,sep = "")
  rowdata<-rowdata%>%bind_rows(read_csv(listfile, col_names = c("EO", "cnt"))%>%mutate("chr" = origin, "kmer"=k))
}
observed<-rowdata%>%filter(EO=="kmer")
df<-rowdata%>%group_by(chr)%>%group_by(kmer)%>%filter(EO!="kmer")%>%summarise("cnt"=mean(cnt))%>%mutate("chr"=origin)

for (j in c(2, 11)){
  hor = j
  k = 15
  listfile = str_c(k, "mer_chr", hor, "HOR_randomcount.csv",sep = "")
  origin = str_c("chr", hor,sep = "")
  rowdata <-read_csv(listfile, col_names = c("EO", "cnt"))%>%mutate("chr" = origin, "kmer"=k)
  rowdata
  
  for (i in 16:21){
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
gp<-ggplot(df)+ geom_line(aes(x=kmer, y=cnt, color=EO))+facet_grid(chr~., scales="free")
print(gp)
ggsave("kmerhit_EO_zoom.png", gp, width=10, height=10)

