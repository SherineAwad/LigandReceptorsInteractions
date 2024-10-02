library(Signac)
library(Seurat)
library(SeuratDisk)
library(GenomeInfoDb)
library(future)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
set.seed(1234)
args <- commandArgs(trailingOnly = TRUE)
mysample = args[1]
myRDS <- paste(mysample, ".rds", sep="")
myh5Seurat <- paste(mysample, ".h5Seurat", sep="") 
myRDS
myObject <- readRDS(myRDS)
myObject[["RNA"]] <- as(object = myObject[["RNA"]], Class = "Assay")
SaveH5Seurat(myObject, filename = myh5Seurat, overwrite= TRUE)
Convert(myh5Seurat, dest = "h5ad", assay="RNA",overwrite = TRUE)



