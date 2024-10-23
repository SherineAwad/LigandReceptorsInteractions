#! /usr/bin/env python
import sys
import argparse
import math
import scanpy as sc
import liana as li
import omnipath as op
import decoupler as dc
import pandas as pd
from liana.method import singlecellsignalr, connectome, cellphonedb, natmi, logfc, cellchat, geometric_mean
from liana.mt import rank_aggregate
from orthologous import  toZebra



def main():
     zfish = toZebra()
     adata = sc.read_h5ad("Injury_snRNA_allClusters_renamed.h5ad")
     print(adata.var_names)
     print(adata.var)
     li.mt.show_methods()

     #new_cellype we have are {'AC', 'MGPCs', 'Rod_p', 'BC', 'RestMG', 'RGC_p', 'RGC', 'Rod', 'AC_p', 'HC', 'ActivateMG', 'BC_p', 'Cone_p', 'Microglia', 'Cone'}
     #we use inplace=False to print to file 
     phonePD = cellphonedb(adata,groupby='new_celltype',resource_name='consensus',resource=zfish,expr_prop=0,  use_raw = False,inplace = False)
     phonePD.to_csv("cellphone.csv", sep='\t', encoding='utf-8', index=False, header=True) 

     cellchatPD = cellchat(adata,groupby='new_celltype',resource_name='consensus',resource=zfish,expr_prop=0, use_raw = False, inplace = False)
     cellchatPD.to_csv("cellchat.csv", sep='\t', encoding='utf-8', index=False, header=True) 
     
     scPD = singlecellsignalr(adata,groupby='new_celltype',resource_name='consensus',resource=zfish,expr_prop=0, use_raw = False,inplace = False)
     scPD.to_csv("signalR.csv", sep='\t', encoding='utf-8', index=False, header=True)
     
     rkm = rank_aggregate(adata,groupby='new_celltype',resource_name='consensus',resource=zfish,expr_prop=0,  aggregate_method="mean", use_raw = False, inplace=False)
     rkm.to_csv("RK_mean.csv",  sep='\t', encoding='utf-8', index=False, header=True)
     
     rkr = rank_aggregate(adata,groupby='new_celltype',resource_name='consensus',resource=zfish,expr_prop=0,  aggregate_method="rra", use_raw = False, inplace=False)
     rkr.to_csv("RK_rra.csv",  sep='\t', encoding='utf-8', index=False, header=True)


if __name__ == '__main__':
    main()



