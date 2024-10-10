
rule all:
    input:
       "cellchatotm1.csv",  
       "cellphoneotm1.csv",  
       "signalRotm1.csv"

rule LR:
   input: 
      h5ad = config['h5ad'] 
   output: 
       "cellchatotm1.csv",  
       "cellphoneotm1.csv",  
       "signalRotm1.csv"
   shell:
      "python LR.py"  
