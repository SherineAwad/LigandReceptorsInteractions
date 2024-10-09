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
from liana.method import singlecellsignalr, connectome, cellphonedb, natmi, logfc, cellchat, geometric_mean


def toZebra():
     resource = li.rs.select_resource('consensus')
     resource.head()
     #we can use other nameID as in [human_entrez_gene	human_ensembl_gene	hgnc_id	human_name	human_symbol	human_chr	human_assert_ids	zebrafish_entrez_gene	zebrafish_ensembl_gene	zfin_id	zebrafish_name	zebrafish_symbol	zebrafish_chr	zebrafish_assert_ids]
     map_df = li.rs.get_hcop_orthologs(url='https://ftp.ebi.ac.uk/pub/databases/genenames/hcop/human_zebrafish_hcop_fifteen_column.txt.gz',
                                   columns=['human_symbol', 'zebrafish_symbol'],
                                   # NOTE: HCOP integrates multiple resource, so we can filter out mappings in at least 3 of them for confidence
                                   min_evidence=3
                                   )

     # rename the columns to source and target, respectively for the original organism and the target organism
     map_df = map_df.rename(columns={'human_symbol':'source', 'zebrafish_symbol':'target'})
     map_df.tail()

     zfish = li.rs.translate_resource(resource,
                                 map_df=map_df,
                                 columns=['ligand', 'receptor'],
                                 replace=True,
                                 # NOTE that we need to define the threshold of redundancies for the mapping
                                 # in this case, we would keep mappings as long as they don't map to more than 2 zebrafish genes
                                 one_to_many=3
                                 )
     return zfish


def main():
     zfish = toZebra()
     with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
         zfish.to_csv('zebraEns_LR.csv', sep='\t', index=False, encoding='utf-8')
     
     adata = sc.read_h5ad("Injury_snRNA_allClusters.h5ad")
     for line in open("mapping.csv"): #if we use Ensembl ID we need mappignEnsmbl.csv 
          if "old" in line: 
             continue
          records = (line.strip()).split(",")
          adata.var_names = adata.var_names.str.replace(records[0], records[1])
     #adata.raw.X
     print(adata.var_names)
     print(adata.var)
     li.mt.show_methods()

     #new_cellype we have are {'AC', 'MGPCs', 'Rod_p', 'BC', 'RestMG', 'RGC_p', 'RGC', 'Rod', 'AC_p', 'HC', 'ActivateMG', 'BC_p', 'Cone_p', 'Microglia', 'Cone'}
     #we use inplace=False to print to file 
     phonePD = cellphonedb(adata,groupby='new_celltype',resource_name='consensus',resource=zfish,expr_prop=0, key_added='cpdb_res', use_raw = False,inplace = False)
     phonePD.to_csv("cellphonedbLR.csv", sep='\t', encoding='utf-8', index=False, header=True) 
     cellchatPD = cellchat(adata,groupby='new_celltype',resource_name='consensus',resource=zfish,expr_prop=0, key_added='cchat_res',use_raw = False, inplace = False)
     cellchatPD.to_csv("cellchatdbLR.csv", sep='\t', encoding='utf-8', index=False, header=True) 
     scPD = singlecellsignalr(adata,groupby='new_celltype',resource_name='consensus',resource=zfish,expr_prop=0, key_added='scdb_res', use_raw = False,inplace = False)
     scPD.to_csv("singlecelldbLR.csv", sep='\t', encoding='utf-8', index=False, header=True)
     
     
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



     rkm = rank_aggregate(adata,groupby='new_celltype',resource_name='consensus',resource=zfish,expr_prop=0,  aggregate_method="mean",key_added ='mean_rank', use_raw = False, inplace=False)
     rkm.to_csv("RK_mean.csv",  sep='\t', encoding='utf-8', index=False, header=True)
     rkr = rank_aggregate(adata,groupby='new_celltype',resource_name='consensus',resource=zfish,expr_prop=0,  aggregate_method="rra",key_added ='rra_rank', use_raw = False, inplace=False)
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
              orderby='specificity_rank',
              orderby_ascending=True,
              figure_size=(25, 15)
             )
     p4.save("aggregate.png")

if __name__ == '__main__':
    main()



