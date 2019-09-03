library(data.table)
library(tidyverse)
library(plyr)
library(reshape2)
library(compare)
setwd("mydna")
getwd()

mydata<-read.csv("C:/Users/SamBr/Documents/pgx/mygenotypes.csv")%>%
  select(-X)%>%
  rename(c('gene' = 'GENE', 'genotype' = 'GENOTYPE'))
  


files <- list.files()
temp <- lapply(files, fread, sep=",")
data <- rbindlist(temp)%>%
  mutate_if(is.character, str_trim)




write_csv (data, "C:/Users/SamBr/Documents/pgx/all_mydna.csv")

xx<-split(data, data$GENE, drop = F)
#names(xx)
#CYP2D6<-xx$CYP2D6



#v<-table(CYP2D6$GENOTYPE)%>%
# as.data.frame()%>%
# rename(c('Var1' = 'GENOTYPE'))%>%
# mutate(per_freq = Freq/150*100)%>%
# merge(CYP2D6, by = 'GENOTYPE')


genotype_freqs<-function(df){
  v<-table(df$GENOTYPE)%>%
  as.data.frame()%>%
  rename(c('Var1' = 'GENOTYPE'))%>%
  mutate(per_freq = Freq/150*100)%>%
  merge(df, by = 'GENOTYPE')}


c<-genotype_freqs(CYP2D6)
r<-lapply(xx, genotype_freqs)

data_F <- rbindlist(r)

data_Freq <- rbindlist(r)%>%
  select(-ID,-`PREDICTED FUNCTION`)%>%
  unique()


frequencies<-merge(mydata, data_Freq, by = c("GENE","GENOTYPE"))

frequencies<- format(frequencies, digits = 2)

write_csv (frequencies, "C:/Users/SamBr/Documents/pgx/myfreq.csv")

kdata<-read.csv("C:/Users/SamBr/Documents/pgx/allkailos3.csv")

kd<-kdata %>% mutate_if(is.factor, as.character)%>%
  mutate_if(is.character, str_trim)


m<-melt(kd, id.vars = "ID", measure.vars = 2:39)%>%
  mutate_all(~replace(., . == "", NA)) %>%
  na.omit()%>%
  rename(c('variable' = 'GENE', 'value'='kailos_call'))

lol<-table(m$kailos_call)%>%
  as.data.frame()%>%
  rename(c('Var1' = 'kailos_call'))%>%
  mutate(k_per_freq = Freq/37*100)%>%
  rename(c('Freq' = 'k_Freq'))%>%
  merge(m, by = 'kailos_call')%>%
  select(-ID)%>%
  unique()%>%
  format(digits = 2)

mykailos<-merge(frequencies, lol, by = "kailos_call")
mykailos<-join(frequencies, lol, by = "kailos_call", match = "all")

POP<-unique(lol)
