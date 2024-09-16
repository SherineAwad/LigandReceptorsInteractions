#! /usr/bin/env python
import sys
import argparse
import math


def readGenes(genefile): 
    newNames = []
    print("oldName,newName")
    for line in open(genefile):
        records = (line.strip()).split("~~")
        if records[1] in newNames: 
            continue 
        else: 
            newNames.append(records[1])
        print(line.strip()+ ","+records[1]) 
    return 

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('genefile')
    args = parser.parse_args()
    genefile = args.genefile  
    readGenes(genefile) 
    
if __name__ == '__main__':
    main()

   

