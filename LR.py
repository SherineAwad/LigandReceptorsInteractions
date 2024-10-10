! /usr/bin/env python
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
     phonePD.to_csv("cellphoneotm1.csv", sep='\t', encoding='utf-8', index=False, header=True) 
     cellchatPD = cellchat(adata,groupby='new_celltype',resource_name='consensus',resource=zfish,expr_prop=0, use_raw = False, inplace = False)
     cellchatPD.to_csv("cellchatotm1.csv", sep='\t', encoding='utf-8', index=False, header=True) 
     scPD = singlecellsignalr(adata,groupby='new_celltype',resource_name='consensus',resource=zfish,expr_prop=0, use_raw = False,inplace = False)
     scPD.to_csv("signalRotm1.csv", sep='\t', encoding='utf-8', index=False, header=True)
     
     
     #we run inplace =True for plots     
     cellphonedb(adata,groupby='new_celltype',resource_name='consensus',resource=zfish,expr_prop=0, key_added='cpdb_res', use_raw = False)
     p1 = li.pl.dotplot(adata = adata,
              colour='lr_means',
              size='cellphone_pvals',
              inverse_size=True, # we inverse sign since we want small p-values to have large sizes
              source_labels=['Microglia'],
              target_labels=['MGPCs','Rod', 'Cone', 'RGC','ActivateMG'],
              figure_size=(20, 10),
              filter_fun=lambda x: x['cellphone_pvals'] <= 0.001,
              top_n=10,
              orderby='cellphone_pvals',
              orderby_ascending=True,
              uns_key='cpdb_res' # uns_key to use, default is 'liana_res'
             )
     p1.save("cellphone.png")

     cellchat(adata,groupby='new_celltype',resource_name='consensus',resource=zfish,expr_prop=0, key_added='cchat_res',use_raw = False) 
     
     p2= li.pl.dotplot(adata = adata,
              colour='lr_probs',
              size='cellchat_pvals',
              inverse_size=True, # we inverse sign since we want small p-values to have large sizes
              source_labels=['Microglia'],
              target_labels=['MGPCs','Rod', 'Cone', 'RGC','ActivateMG'],
              figure_size=(25, 10),
              filter_fun=lambda x: x['cellchat_pvals'] <= 0.001,
              top_n=10,
              orderby='cellchat_pvals',
              orderby_ascending=True,
              uns_key='cchat_res' # uns_key to use, default is 'liana_res'
             )
     p2.save("cellchat.png")
     
     singlecellsignalr(adata,groupby='new_celltype',resource_name='consensus',resource=zfish,expr_prop=0, key_added='scdb_res', use_raw = False) 
     
     p3= li.pl.dotplot(adata = adata,
              colour='lrscore',
              size='lrscore',
              inverse_size=True, # we inverse sign since we want small p-values to have large sizes
              source_labels=['Microglia'],
              target_labels=['MGPCs','Rod', 'Cone', 'RGC','ActivateMG'],
              figure_size=(25, 10),
              top_n=10,
              orderby='lrscore',
              orderby_ascending=True,
              uns_key='scdb_res' # uns_key to use, default is 'liana_res'
             )
     p3.save("signalR.png")



     rkm = rank_aggregate(adata,groupby='new_celltype',resource_name='consensus',resource=zfish,expr_prop=0,  aggregate_method="mean", use_raw = False, inplace=False)
     rkm.to_csv("RK_mean.csv",  sep='\t', encoding='utf-8', index=False, header=True)
     rkr = rank_aggregate(adata,groupby='new_celltype',resource_name='consensus',resource=zfish,expr_prop=0,  aggregate_method="rra", use_raw = False, inplace=False)
     rkr.to_csv("RK_rra.csv",  sep='\t', encoding='utf-8', index=False, header=True)

     rank_aggregate(adata,groupby='new_celltype',resource_name='consensus',resource=zfish,expr_prop=0,  aggregate_method="rra",key_added ='rra_rank', use_raw = False) 
     p4 = li.pl.dotplot(adata = adata,
              colour='specificity_rank',
              size='spec_weight',
              inverse_size=True,
              inverse_colour=True,
              source_labels=['ActivateMG'],
              target_labels=['MGPCs','Rod', 'Cone', 'RGC','ActivateMG'], 
              top_n=10,
              uns_key='rra',
              orderby='specificity_rank',
              orderby_ascending=True,
              figure_size=(25, 15)
             )
     p4.save("aggregate.png")

if __name__ == '__main__':
    main()



