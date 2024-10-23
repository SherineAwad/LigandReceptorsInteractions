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
                                 ## in this case, we would keep mappings as long as they don't map to more than 2 zebrafish genes
                                 #one_to_many=3
                                 #will keep that map only to one gene 
                                 one_to_many=5
                                 )
     return zfish


def main():
     zfish = toZebra()
     with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
         zfish.to_csv('zebrafish_LR.csv', sep='\t', index=False, encoding='utf-8')
     
if __name__ == '__main__':
    main()



