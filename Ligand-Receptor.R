library(tidyverse)
library(magrittr)
library(liana)
library(OmnipathR)
library(readr)

#Require conda activate archr 
show_resources()

show_methods()

liana_path <- system.file(package = "liana")

args <- commandArgs(trailingOnly = TRUE)
mysample <- args[1]

myRDS <- paste(mysample, ".rds", sep="")
myRDS

mydata <- readRDS(myRDS)


zebra_human <- homologene_download(
  target = 9606L,
  source = 7955L,
  id_type = "genesymbol",
  hgroup_size = TRUE
)

zebra_human <- rename(zebra_human,  source_genesymbol =genesymbol_source, target_genesymbol=genesymbol_target )
head(zebra_human)

mydata %>% dplyr::glimpse() 
#Could use assay =RNA
liana_test <- liana_wrap(mydata, assay='RNA', resource ='custom', external_resource=zebra_human)
liana_test %>% dplyr::glimpse()
liana_test <- liana_test %>% liana_aggregate()
dplyr::glimpse(liana_test)

figure_name <- mysample
figure_name <- paste(figure_name,"LR1.pdf", sep="")
pdf(file =figure_name, height=20, width=40)
liana_test %>%
  liana_dotplot(source_groups = c("MG"),
                target_groups = c("PreAC", "AC", "PreBC", "BC","Pre Rod", "Rod", "MG"), 
                ntop = 10)

dev.off() 

figure_name <- mysample
figure_name <- paste(figure_name,"LR2.pdf", sep="")
pdf(file =figure_name, height=20, width=40)
liana_test %>%
  liana_dotplot(source_groups = c("MG"),
                target_groups = c("RGC", "Pre Cone" ,"Cone", "Stage1","Stage2" , "Act_MG","MG"),
                ntop = 10)

dev.off() 

