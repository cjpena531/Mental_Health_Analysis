import os 
import json
import sys
sys.path.append('src')
from etl import *
import subprocess


def main():
    
    args = sys.argv
    
    if (len(args) > 1):
        if (args[1] == "test"):
            print("testing...")
            dictionary = json.load(open("config/test_parameters.json"))
        elif (args[1] == "clean"):
            print("cleaning...")
            clean()
            print("cleaned!")
            return
    else:
        dictionary = json.load(open("config/parameters.json"))
    
    create_folders()
    subset = get_subset(dictionary['read_1'],dictionary['read_2'],dictionary['subset_size'])
    #fastqc(dictionary,subset)
    #cutadapt(dictionary,subset)
    #kallisto(dictionary,subset)
    gene_counts(dictionary,subset)
    pca()
    
    
if __name__ == '__main__':
    main()