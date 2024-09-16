
rule all:
    input:
       "InjuryZebraLR1.pdf" ,
       "InjuryZebraLR2.pdf"

rule split: 
    input: 
        "genes.txt" 
    output: 
         "mapping.csv" 
    shell: 
        "python split.py {input} > {output}"

rule rename: 
    input: 
       "mapping.csv", 
    params: 
       myObject =config['SeuratObject']  
    output: 
        "InjuryZebra.rds"    
    shell: 
       "Rscript renameFeatures.R  {params}"

rule LR:
   input: 
      "InjuryZebra.rds"
   output: 
       "InjuryZebraLR1.pdf" ,
       "InjuryZebraLR2.pdf"
   shell: 
        "Rscript Ligand-Receptor.R  {input}" 
