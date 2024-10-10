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


def rename(mapping_file, h5ad_file, newh5ad_file):
    adata = sc.read_h5ad(h5ad_file)
     for line in open(mapping_file): #mapping.csv or if we use Ensembl ID we need mappignEnsmbl.csv 
          if "old" in line:
             continue
          records = (line.strip()).split(",")
          adata.var_names = adata.var_names.str.replace(records[0], records[1])
     adata.write_h5ad(newh5ad_file)
    return 


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('mapping_file')
    parser.add_argument('h5ad_file')
    parser.add_argumen('newh5ad_file') 
    args = parser.parse_args()
    mapping_file = args.mapping_file
    h5ad_file = args.h5ad_file 
    newh5ad_file = args.newh5ad_file 
    rename(mapping_file, h5ad_file, newh5ad_file)

if __name__ == '__main__':
    main()



