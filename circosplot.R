library(tidyverse)
# install.packages("circlize")
library(circlize)
setwd("~/git/k_mers/")

(df<-read_csv("HOR1_19mer_chr1HOR_hit_test.gff", comment="#", col_names = c("chr", "source", "anno","start", "end","comment", "strand", "kmer"))
  %>%select(chr, start, end, strand, kmer))

# df%>%arrange(kmer)%>%View
df<-df%>%left_join(read_csv("monomer_length.csv"), by="chr")%>%
group_by(kmer)%>%filter(n()>1)%>%ungroup()%>%arrange(kmer)

df<- df%>%filter(start<m1)%>%mutate("monomer"="m1")%>%mutate("linkpoint"=start)%>%
  bind_rows(df%>%filter(start<m2 & start>m1)%>%mutate("monomer"="m2")%>%mutate("linkpoint"=start-m1))%>%
  bind_rows(df%>%filter(start<m3 & start>m2)%>%mutate("monomer"="m3")%>%mutate("linkpoint"=start-m2))%>%
  bind_rows(df%>%filter(start<m4 & start>m3)%>%mutate("monomer"="m4")%>%mutate("linkpoint"=start-m3))%>%
  select(monomer, start, end, linkpoint, kmer)%>%arrange(kmer)

f<-read_csv("monomer_length.csv")
f<-f%>%gather(key="monomer", value="mend", -chr)%>%mutate("mstart"=lag(mend, default=1))
# %>%mutate(m2_x=m2-m1)%>%mutate(m3_x=m3-m2)%>%mutate(m4_x=m4-m3)

(df_link<-df%>%select(monomer, start, kmer)%>%spread(key=monomer, value=start)%>%arrange(kmer))
(df_name<-df%>%select(monomer, kmer)%>%spread(key=monomer, value=monomer)%>%arrange(kmer))
  
(df_names<-df_name%>%filter(is.na(m1)==FALSE&is.na(m3)==FALSE)%>%select(kmer, m1, m3)%>%rename(x1=m1, x2=m3)%>%
    bind_rows(df_name%>%filter(is.na(m1)==FALSE&is.na(m4)==FALSE)%>%select(kmer, m1, m4)%>%rename(x1=m1, x2=m4))%>%
    bind_rows(df_name%>%filter(is.na(m2)==FALSE&is.na(m3)==FALSE)%>%select(kmer, m2, m3)%>%rename(x1=m2, x2=m3))%>%
    bind_rows(df_name%>%filter(is.na(m2)==FALSE&is.na(m4)==FALSE)%>%select(kmer, m2, m4)%>%rename(x1=m2, x2=m4))%>%
    bind_rows(df_name%>%filter(is.na(m3)==FALSE&is.na(m4)==FALSE)%>%select(kmer, m3, m4)%>%rename(x1=m3, x2=m4)))

(df_m1_m2<-df_link%>%filter(is.na(m1)==FALSE&is.na(m2)==FALSE)%>%select(kmer, m1, m2))
(df_m1_m3<-df_link%>%filter(is.na(m1)==FALSE&is.na(m3)==FALSE)%>%select(kmer, m1, m3))
(df_m1_m4<-df_link%>%filter(is.na(m1)==FALSE&is.na(m4)==FALSE)%>%select(kmer, m1, m4))
(df_m2_m3<-df_link%>%filter(is.na(m2)==FALSE&is.na(m3)==FALSE)%>%select(kmer, m2, m3))
(df_m2_m4<-df_link%>%filter(is.na(m2)==FALSE&is.na(m4)==FALSE)%>%select(kmer, m2, m4))
(df_m3_m4<-df_link%>%filter(is.na(m3)==FALSE&is.na(m4)==FALSE)%>%select(kmer, m3, m4))



circos.clear()
circos.par(cell.padding=c(0,0,0,0), track.margin=c(0,0.15), start.degree = 90, gap.degree =4)
circos.initialize(factors = df$monomer, xlim = cbind(f$mstart-1, f$mend))
circos.trackPlotRegion(ylim = c(0, 1), factors = f$monomer, track.height=0.1,
                       panel.fun = function(x, y) {
                         name = get.cell.meta.data("sector.index")
                         i = get.cell.meta.data("sector.numeric.index")
                         xlim = get.cell.meta.data("xlim")
                         ylim = get.cell.meta.data("ylim")
                         circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2], ytop=ylim[2])
                         circos.text(x=mean(xlim), y=1.7, labels=name, cex=0.6)})

for(k in 1:nrow(df_m1_m3)){
  #plot link
  circos.link("m1", point1=c(df_m1_m3$m1[k], df_m1_m3$m1[k]+19),
              "m3", point2=c(df_m1_m3$m3[k], df_m1_m3$m3[k]+19))
}

for(k in 1:nrow(df_m2_m4)){
  #plot link
  circos.link("m1", point1=c(df_m2_m4$m2[k], df_m2_m4$m2[k]+19),
              "m3", point2=c(df_m2_m4$m4[k], df_m2_m4$m4[k]+19))
}






(df<-read_csv("HOR1_10mer_chr1HOR_hit_test.gff", comment="#", col_names = c("chr", "source", "anno","start", "end","comment", "strand", "kmer"))
  %>%select(chr, start, end, strand, kmer))

# df%>%arrange(kmer)%>%View
df<-df%>%left_join(read_csv("monomer_length.csv"), by="chr")%>%
  group_by(kmer)%>%filter(n()>1)%>%ungroup()%>%arrange(kmer)

df<- df%>%filter(start<m1)%>%mutate("monomer"="m1")%>%mutate("linkpoint"=start)%>%
  bind_rows(df%>%filter(start<m2 & start>m1)%>%mutate("monomer"="m2")%>%mutate("linkpoint"=start-m1))%>%
  bind_rows(df%>%filter(start<m3 & start>m2)%>%mutate("monomer"="m3")%>%mutate("linkpoint"=start-m2))%>%
  bind_rows(df%>%filter(start<m4 & start>m3)%>%mutate("monomer"="m4")%>%mutate("linkpoint"=start-m3))%>%
  select(monomer, start, end, linkpoint, kmer)%>%arrange(kmer)

f<-read_csv("monomer_length.csv")
f<-f%>%gather(key="monomer", value="mend", -chr)%>%mutate("mstart"=lag(mend, default=1))
# %>%mutate(m2_x=m2-m1)%>%mutate(m3_x=m3-m2)%>%mutate(m4_x=m4-m3)

(df_link<-df%>%select(monomer, start, kmer)%>%spread(key=monomer, value=start)%>%arrange(kmer))
(df_name<-df%>%select(monomer, kmer)%>%spread(key=monomer, value=monomer)%>%arrange(kmer))

(df_names<-df_name%>%filter(is.na(m1)==FALSE&is.na(m3)==FALSE)%>%select(kmer, m1, m3)%>%rename(x1=m1, x2=m3)%>%
    bind_rows(df_name%>%filter(is.na(m1)==FALSE&is.na(m4)==FALSE)%>%select(kmer, m1, m4)%>%rename(x1=m1, x2=m4))%>%
    bind_rows(df_name%>%filter(is.na(m2)==FALSE&is.na(m3)==FALSE)%>%select(kmer, m2, m3)%>%rename(x1=m2, x2=m3))%>%
    bind_rows(df_name%>%filter(is.na(m2)==FALSE&is.na(m4)==FALSE)%>%select(kmer, m2, m4)%>%rename(x1=m2, x2=m4))%>%
    bind_rows(df_name%>%filter(is.na(m3)==FALSE&is.na(m4)==FALSE)%>%select(kmer, m3, m4)%>%rename(x1=m3, x2=m4)))

(df_m1_m2<-df_link%>%filter(is.na(m1)==FALSE&is.na(m2)==FALSE)%>%select(kmer, m1, m2))
(df_m1_m3<-df_link%>%filter(is.na(m1)==FALSE&is.na(m3)==FALSE)%>%select(kmer, m1, m3))
(df_m1_m4<-df_link%>%filter(is.na(m1)==FALSE&is.na(m4)==FALSE)%>%select(kmer, m1, m4))
(df_m2_m3<-df_link%>%filter(is.na(m2)==FALSE&is.na(m3)==FALSE)%>%select(kmer, m2, m3))
(df_m2_m4<-df_link%>%filter(is.na(m2)==FALSE&is.na(m4)==FALSE)%>%select(kmer, m2, m4))
(df_m3_m4<-df_link%>%filter(is.na(m3)==FALSE&is.na(m4)==FALSE)%>%select(kmer, m3, m4))

df_m1_m2<-df_m1_m2%>%arrange(m1)%>%distinct(m1, m2, .keep_all=TRUE)%>%
  mutate("m1_end"=m1+10)%>%
  mutate("m2_end"=m2+10)%>%
  mutate("startcheck"=((m1-m1%>%lag())!=1))%>%
  mutate("endcheck"=((m1%>%lead()-m1)!=1))
df_m1_m2[is.na(df_m1_m2)] <- TRUE

df_m1_m2_start<-df_m1_m2%>%filter(startcheck==TRUE)%>%select(m1, m2)
df_m1_m2_end<-df_m1_m2%>%filter(endcheck==TRUE)%>%select(m1_end, m2_end)

df_m1_m3<-df_m1_m3%>%arrange(m1)%>%distinct(m1, m3, .keep_all=TRUE)%>%
  mutate("m1_end"=m1+10)%>%
  mutate("m3_end"=m3+10)%>%
  mutate("startcheck"=((m1-m1%>%lag())!=1))%>%
  mutate("endcheck"=((m1%>%lead()-m1)!=1))
df_m1_m3[is.na(df_m1_m3)] <- TRUE

df_m1_m3_start<-df_m1_m3%>%filter(startcheck==TRUE)%>%select(m1, m3)
df_m1_m3_end<-df_m1_m3%>%filter(endcheck==TRUE)%>%select(m1_end, m3_end)


df_m1_m4<-df_m1_m4%>%arrange(m1)%>%distinct(m1, m4, .keep_all=TRUE)%>%
  mutate("m1_end"=m1+10)%>%
  mutate("m4_end"=m4+10)%>%
  mutate("startcheck"=((m1-m1%>%lag())!=1))%>%
  mutate("endcheck"=((m1%>%lead()-m1)!=1))
df_m1_m4[is.na(df_m1_m4)] <- TRUE

df_m1_m4_start<-df_m1_m4%>%filter(startcheck==TRUE)%>%select(m1, m4)
df_m1_m4_end<-df_m1_m4%>%filter(endcheck==TRUE)%>%select(m1_end, m4_end)


df_m2_m3<-df_m2_m3%>%arrange(m2)%>%distinct(m2, m3, .keep_all=TRUE)%>%
  mutate("m2_end"=m2+10)%>%
  mutate("m3_end"=m3+10)%>%
  mutate("startcheck"=((m2-m2%>%lag())!=1))%>%
  mutate("endcheck"=((m2%>%lead()-m2)!=1))
df_m2_m3[is.na(df_m2_m3)] <- TRUE

df_m2_m3_start<-df_m2_m3%>%filter(startcheck==TRUE)%>%select(m2, m3)
df_m2_m3_end<-df_m2_m3%>%filter(endcheck==TRUE)%>%select(m2_end, m3_end)


df_m2_m4<-df_m2_m4%>%arrange(m2)%>%distinct(m2, m4, .keep_all=TRUE)%>%
  mutate("m2_end"=m2+10)%>%
  mutate("m4_end"=m4+10)%>%
  mutate("startcheck"=((m2-m2%>%lag())!=1))%>%
  mutate("endcheck"=((m2%>%lead()-m2)!=1))
df_m2_m4[is.na(df_m2_m4)] <- TRUE

df_m2_m4_start<-df_m2_m4%>%filter(startcheck==TRUE)%>%select(m2, m4)
df_m2_m4_end<-df_m2_m4%>%filter(endcheck==TRUE)%>%select(m2_end, m4_end)


df_m3_m4<-df_m3_m4%>%arrange(m3)%>%distinct(m3, m4, .keep_all=TRUE)%>%
  mutate("m3_end"=m3+10)%>%
  mutate("m4_end"=m4+10)%>%
  mutate("startcheck"=((m3-m3%>%lag())!=1))%>%
  mutate("endcheck"=((m3%>%lead()-m3)!=1))
df_m3_m4[is.na(df_m3_m4)] <- TRUE

df_m3_m4_start<-df_m3_m4%>%filter(startcheck==TRUE)%>%select(m3, m4)
df_m3_m4_end<-df_m3_m4%>%filter(endcheck==TRUE)%>%select(m3_end, m4_end)


circos.clear()
circos.par(cell.padding=c(0,0,0,0), track.margin=c(0,0.15), start.degree = 90, gap.degree =4)
circos.initialize(factors = df$monomer, xlim = cbind(f$mstart-1, f$mend))
circos.trackPlotRegion(ylim = c(0, 1), factors = f$monomer, track.height=0.1,
                       panel.fun = function(x, y) {
                         name = get.cell.meta.data("sector.index")
                         i = get.cell.meta.data("sector.numeric.index")
                         xlim = get.cell.meta.data("xlim")
                         ylim = get.cell.meta.data("ylim")
                         circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2], ytop=ylim[2])
                         circos.text(x=mean(xlim), y=1.7, labels=name, cex=0.6)})

for(k in 1:nrow(df_m1_m2_start)){
  #plot link
  circos.link("m1", point1=c(df_m1_m2_start$m1[k], df_m1_m2_end$m1_end[k]),
              "m2", point2=c(df_m1_m2_start$m2[k], df_m1_m2_end$m2_end[k]), col="#FF000080")
}

for(k in 1:nrow(df_m1_m3_start)){
  #plot link
  circos.link("m1", point1=c(df_m1_m3_start$m1[k], df_m1_m3_end$m1_end[k]),
              "m3", point2=c(df_m1_m3_start$m3[k], df_m1_m3_end$m3_end[k]), col="#FF000080")
}

for(k in 1:nrow(df_m1_m4_start)){
  #plot link
  circos.link("m1", point1=c(df_m1_m4_start$m1[k], df_m1_m4_end$m1_end[k]),
              "m4", point2=c(df_m1_m4_start$m4[k], df_m1_m4_end$m4_end[k]), col="#FF000080")
}

for(k in 1:nrow(df_m2_m3_start)){
  #plot link
  circos.link("m2", point1=c(df_m2_m3_start$m2[k], df_m2_m3_end$m2_end[k]),
              "m3", point2=c(df_m2_m3_start$m3[k], df_m2_m3_end$m3_end[k]), col="#FF000080")
}


for(k in 1:nrow(df_m2_m4)){
  #plot link
  circos.link("m2", point1=c(df_m2_m4_start$m2[k], df_m2_m4_end$m2_end[k]),
              "m4", point2=c(df_m2_m4_start$m4[k], df_m2_m4_end$m4_end[k]), col="#FF000080")
}

for(k in 1:nrow(df_m3_m4)){
  #plot link
  circos.link("m3", point1=c(df_m3_m4_start$m3[k], df_m3_m4_end$m3_end[k]),
              "m4", point2=c(df_m3_m4_start$m4[k], df_m3_m4_end$m4_end[k]), col="#FF000080")
}



