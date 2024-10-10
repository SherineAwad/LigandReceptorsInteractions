#This makefile is to have some changes for the data before running the snakemake, 
#this is only a special case for this dataset, so we have a seprate snakefile, and this makefile is to make the extra preprocessing for this specific case 


#convert our seurat oobject to h5ad 
Injury_snRNA_allClusters.h5ad:
	Rscript seuratoh5ad.R Injury_snRNA_allClusters 

#rename our h5ad object features
Injury_snRNA_allClusters_renamed.h5ad: 
	python renameFeatures.py mapping.csv Injury_snRNA_allClusters.h5ad Injury_snRNA_allClusters_renamed.h5ad
