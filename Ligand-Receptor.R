library(tidyverse)
library(magrittr)
library(liana)

show_resources()

show_methods()

liana_path <- system.file(package = "liana")

args <- commandArgs(trailingOnly = TRUE)
mysamples <- args[1]

myRDS <- paste(mysample, ".rds", sep="")
myRDS

mydata <- readRDS(myRDS)

mydata %>% dplyr::glimpse() 
#Could use assay =RNA 
liana_test <- liana_wrap(mydata, assay='SCT')
liana_test %>% dplyr::glimpse()
liana_test <- liana_test %>% liana_aggregate()
dplyr::glimpse(liana_test)

figure_name <- mysamples
figure_name <- paste(figure_name,"LR.pdf", sep="")
pdf(file =figure_name, height=20, width=40)
liana_test %>%
  liana_dotplot(source_groups = c("Mullerglia"),
                target_groups = c("HC","Cone","Rod","Rod BC","Cone BC","GABAergic AC","Glycinergic/Starbust AC","RGC", "Mullerglia"),
                ntop = 30)

dev.off() 


figure_name <- mysamples
figure_name <- paste(figure_name,"LRHeatmap.pdf", sep="")
pdf(file =figure_name, width=12)
#only keep interactions concordant between methods
#note that these pvals are already corrected
liana_trunc <- liana_test %>% filter(aggregate_rank <= 0.01)
heat_freq(liana_trunc)
dev.off() 


figure_name <- mysamples
figure_name <- paste(figure_name,"Chord.pdf", sep="")
pdf(file =figure_name, width=12)
p <- chord_freq(liana_trunc,
                source_groups = c(""Mullerglia"),
                target_groups = c("Rod", "Cone", "Mullerglia"))

dev.off() 

