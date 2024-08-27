library(tidyverse)
library(magrittr)
library(liana)

show_resources()

show_methods()

liana_path <- system.file(package = "liana")

mydata <- readRDS("merged3samples.rds")

mydata %>% dplyr::glimpse() 
#liana_test <- liana_wrap(mydata)
liana_test <- liana_wrap(mydata, assay='SCT')
liana_test %>% dplyr::glimpse()
liana_test <- liana_test %>% liana_aggregate()
dplyr::glimpse(liana_test)

figure_name <- "merged3samples"
figure_name <- paste(figure_name,"LR.pdf", sep="")
pdf(file =figure_name, width=12)
liana_test %>%
  liana_dotplot(source_groups = c("Rod"),
                target_groups = c("Cone", "Mullerglia", "Rod"),
                ntop = 20)

dev.off() 
figure_name <- "merged3samples"
figure_name <- paste(figure_name,"LRHeatmap.pdf", sep="")
pdf(file =figure_name, width=12)
liana_trunc <- liana_test %>%
   # only keep interactions concordant between methods
  filter(aggregate_rank <= 0.01) # note that these pvals are already corrected

heat_freq(liana_trunc)
dev.off() 



