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


def toZebra():
     resource = li.rs.select_resource('consensus')
     resource.head()
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
         zfish.to_csv('zebra_LR.csv', sep='\t', index=False, encoding='utf-8')
     
     adata = sc.read_h5ad("Injury_snRNA_allClusters.h5ad")

     adata.raw.X

     li.mt.show_methods()

     from liana.mt import rank_aggregate

     from liana.method import singlecellsignalr, connectome, cellphonedb, natmi, logfc, cellchat, geometric_mean

     cellphonedb(adata,groupby='sample',resource_=zfish,expr_prop=0.1,verbose=True, key_added='cpdb_res')



if __name__ == '__main__':
    main()



