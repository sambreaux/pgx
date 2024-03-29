library(data.table)
library(tidyverse)
library(plyr)
library(reshape2)
library(compare)
library(rowr)
library(genetics)

setwd("C:/Users/SamBr/Documents/pgx/mydna")
getwd()



snpdat<-read.csv("C:/Users/SamBr/Documents/pgx/mydna snps - Sheet1.csv")
mydata<-read.csv("C:/Users/SamBr/Documents/pgx/mygenotypes.csv")%>%
  rename(c('gene' = 'GENE', 'genotype' = 'GENOTYPE'))
  


files <- list.files()
temp <- lapply(files, fread, sep=",")
data <- rbindlist(temp)%>%
  mutate_if(is.character, str_trim)




write_csv (data, "C:/Users/SamBr/Documents/pgx/results/tidy_mydna.csv")

xx<-split(data, data$GENE, drop = F)
#names(xx)
#CYP2D6<-xx$CYP2D6   



#v<-table(CYP2D6$GENOTYPE)%>%
# as.data.frame()%>%
# rename(c('Var1' = 'GENOTYPE'))%>%
# mutate(per_freq = Freq/150*100)%>%
# merge(CYP2D6, by = 'GENOTYPE')




###calculate HWB equlibrium for CYP2D6 genes
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

write_csv (frequencies, "C:/Users/SamBr/Documents/pgx/results/myDNA_freq.csv")

kdata<-read.csv("C:/Users/SamBr/Documents/pgx/allkailos3.csv")

kd<-kdata %>% mutate_if(is.factor, as.character)%>%
  mutate_if(is.character, str_trim)


m<-melt(kd, id.vars = "ID", measure.vars = 2:39)%>%
  mutate_all(~replace(., . == "", NA)) %>%
  na.omit()%>%
  rename(c('variable' = 'GENE', 'value'='kailos_call'))

write.csv(m, "C:/Users/SamBr/Documents/pgx/results/tidykailos.csv")

kg<-split(m, m$ID, drop =FALSE)

lol<-table(m$kailos_call)%>%
  as.data.frame()%>%
  rename(c('Var1' = 'kailos_call'))%>%
  mutate(k_per_freq = Freq/37*100)%>%
  rename(c('Freq' = 'k_Freq'))%>%
  merge(m, by = 'kailos_call')%>%
  select(-ID)%>%
  unique()%>%
  format(digits = 2)

write.csv(diffrentsamples, "C:/Users/SamBr/Documents/pgx/results/kailos_freq.csv")


ff<-split(m, m$GENE, drop =FALSE)
md<-mydata$GENE
mkg<-ff[c('CYP2D6', 'CYP2C19', 'CYP3A4', 'CYP3A5', 'CYP2C9')]
remove_rs<-function(x){x[!grepl("rs", x$kailos_call),]}
rs_rem_genes<-lapply(mkg, remove_rs)

k_call_only<-function(x){filter(x, grepl(paste(toMatch, collapse="|"), kailos_call))}

toMatch<-c('rs9923231', 'rs762551', 'rs4149056','CYP2D6', 'CYP2C19', 'CYP3A4', 'CYP3A5', 'CYP2C9')
k_calls_only<-filter(m, grepl(paste(toMatch, collapse="|"), kailos_call))


kailosco_freq <-table(k_calls_only$kailos_call)%>%
  as.data.frame()%>%
  rename(c('Var1' = 'kailos_call'))%>%
  mutate(k_per_freq = Freq/37*100)%>%
  rename(c('Freq' = 'k_Freq'))%>%
  merge(m, by = 'kailos_call')%>%
  select(-ID)%>%
  unique()%>%
  format(digits = 2) 


mydata_geno<-merge(mydata, data, by = c("GENE","GENOTYPE"))
my_v_k<-merge( k_calls_only,mydata_geno, by = c("GENE","ID"), all = T) 



#t60-163 and t60-169 failed mydna genotyping
# think L65-294 and L65-295 were swapped L65-293 =L65-296?
mmm<-my_v_k[!is.na(my_v_k$kailos_call.x),]%>%
  select( -GENOTYPE, -population.level, -phenotype, -'PREDICTED FUNCTION')%>%
  rename(c("kailos_call.x" = 'kailos', 'kailos_call.y' = "myDNA"))

  
mmm$compare <-ifelse(mmm$kailos == mmm$myDNA,  "x",ifelse(mmm$kailos != mmm$myDNA, "different", "a"))
diffrentsamples<-filter(mmm, compare == "different")

write.csv(diffrentsamples, "C:/Users/SamBr/Documents/pgx/results/different_genotypes.csv")
write.csv(mmm, "C:/Users/SamBr/Documents/pgx/results/my_vs_kaliosgenotypes.csv")

l294 <- filter(mmm, ID == "L65-294")
l295<- filter(mmm, ID == "L65-295")

problemchildren<- cbind(l294,l295)



mykailos<-merge(frequencies, kailosco_freq , by = "kailos_call", all = T)
mykailos<-join(frequencies, lol, by = "kailos_call", match = "all")
write.csv(mykailos, "C:/Users/SamBr/Documents/pgx/results/my_vs_kailos_freq_pop.csv")

myunique<-unique(data$ID)
kailosunique <-unique(m$ID)
nnn<-cbind.fill(myunique, kailosunique, fill ="na")


drugs<-read.csv("C:/Users/SamBr/Documents/pgx/drug-gene-md.csv")%>%
  select(-'CPIC.Publications..PMID.')%>%
  rename(c('Gene' = 'GENE'))

mydrugs<-table(data$`PREDICTED FUNCTION`)%>%
  as.data.frame()%>%
  rename(c('Var1' = 'PREDICTED FUNCTION'))%>%
  mutate(per_fun_freq = Freq/150*100)%>%
  rename(c('Freq'= 'function_freq'))%>%
  merge(data, by = 'PREDICTED FUNCTION')%>%
  select(-ID, -GENOTYPE)%>%
  unique()
  
function_freq_drug<-merge(mydrugs, drugs, by ="GENE")
write.csv(function_freq_drug, "C:/Users/SamBr/Documents/pgx/results/myFUNCTION_drug_freq.csv")

pheno_data<-mydata_geno%>%
  select(phenotype, 'PREDICTED FUNCTION')



fun_pheno<-merge(function_freq_drug, pheno_data, by = 'PREDICTED FUNCTION' )%>%
  unique()
fun_pheno<-fun_pheno %>% mutate_if(is.factor, as.character)%>%
  mutate_if(is.character, str_trim)

cast_pheno<-dcast(fun_pheno, Drug+GENE~phenotype, fun.aggregate = sum, value.var = "function_freq")



#remove/replace warfin columns
  #Normal warfarin sensitivity = Normal metaboliser
  #High warfarin sensitivity = Poor metaboliser
  #Increased warfarin sensitivity = Reduced metaboliser

#compare gene to RT data
  #added -rt to kailos 

setwd("C:/Users/samBr/Documents/pgx/myDNA_drugs_CLEAN2")
files2 <- list.files()
temp2 <- lapply(files2, fread, sep=",")
data2 <- as.data.frame(rbindlist(temp2, use.names = TRUE))%>%
  mutate_if(is.character, str_trim)
sites<-separate(data2,ID, into = c("site", "xxx"), sep = "-", remove = F )%>%
  select(-xxx)
write.csv(sites, "C:/Users/SamBr/documents/pgx/results/allmydrugdata.csv")

conc_freq<-as.data.frame(table(sites$site, sites$CONSIDERATION, sites$MEDICATION))%>%
  rename(c('Var1'= 'site', 'Var2' = 'CONSIDERATION', 'Var3'= 'MEDICATION'))%>%
  merge(sites, by = c('site', 'CONSIDERATION', 'MEDICATION'))%>%
  rename(c('GENE(S)\rINVOLVED'= 'GENEs'))

cast_site<-dcast(conc_freq, MEDICATION+GENEs+CONSIDERATION~site, fun.aggregate = mean, value.var = "Freq")%>%
  rename(c('BCP'= 'BCP (n=8)', "D59"  ='D59 (n=1)', "F31" = 'F31 (n=2)' ,   "J02" = "J02 (n=8)", "K74" = 'k74 (n=6)', "K76" = "k76 (n=10)" ,"L65" = "L65 (n=20)",
           "N24"= "N24 (n=20)" , "N43"= "N43 (n=9)", "N83" = "N83 (n=16)", "S73"= "S73 (n=5)", "S81"= "S81 (n=4)", "S94" = "S94 (n=3)", "T37"= "T37 (n=9)", "T60" = "T60 (n=9)", 
           "W19" = "W19 (n=10)", "X73" = "X73 (n=10)"))

cast_site[is.na(cast_site)] <- 0
write.csv(cast_site, "C:/Users/SamBr/documents/pgx/results/site_consideration.csv")

g<-data2$ID%>%
  unique()%>%
  view()


#calc hardy winberg
# FIX KAILOS  l65 
#find misisng smaples in thsoe 2 datasets
#drug class
