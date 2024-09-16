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
liana_test <- liana_test %>% liana_aggregate(resource=resource)
dplyr::glimpse(liana_test)

figure_name <- mysamples
figure_name <- paste(figure_name,"LR.pdf", sep="")
pdf(file =figure_name, height=20, width=40)
liana_test %>%
  liana_dotplot(source_groups = c("Act_MG"),
                target_groups = c("AC", "RGC", "PreAC", "BC", "Act_MG",),
                ntop = 30)

dev.off() 

RGC    BC     clstn2         CLSTN2               0.00000126      104.
 2 AC     Act_MG cpne8          CPNE8                0.0000143       107.
 3 Act_MG AC     cacna2d1a      CACNA2D1             0.0000196       104.
 4 RGC    RGC    nrxn2a         NRXN2                0.0000657       142.
 5 AC     Stage1 cpne8          CPNE8                0.0000954       121.
 6 PreAC  Act_MG cpne8          CPNE8                0.000133        119.
 7 BC     BC     clstn2         CLSTN2               0.000155        138.
 8 RGC    MG     adam12         ADAM12               0.000179        144.
 9 PreAC  AC     cacna2d3       CACNA2D3             0.000206        165.
10 Act_MG PreAC



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

