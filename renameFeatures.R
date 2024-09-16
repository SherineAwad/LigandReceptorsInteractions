library(readr)


args <- commandArgs(trailingOnly = TRUE)
mysample <- args[1]

myRDS <- paste(mysample, ".rds", sep="")
myRDS

mydata <- readRDS(myRDS)



mapping <- read_csv("mapping.csv")
feature_name_mapping <- setNames(mapping$newName, mapping$oldName)
rna_assay <- mydata[["RNA"]]

# Replace feature names in the expression matrix
new_feature_names <- sapply(rownames(rna_assay), function(x) {
    if (x %in% names(feature_name_mapping)) {
        feature_name_mapping[x]
    } else {
        x
    }
})

# Update the feature names in the assay
RNA <- mydata@assays$RNA
RNA@counts@Dimnames[[1]] = new_feature_names
RNA@data@Dimnames[[1]] = new_feature_names
mydata@assays$RNA <- RNA

myRDS <- "InjuryZebra.rds" 

saveRDS(mydata, file=myRDS)

