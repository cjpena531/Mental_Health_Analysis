import os 
import json
import sys
sys.path.insert(0,sys.path[0] +'/src')
from etl import *

def main():
    
    args = sys.argv
    print(args)
    
    if (len(args) > 1):
        if (args[1] == "test"):
            print("yes")
            dictionary = json.load(open("config/test_parameters.json"))
    else:
        dictionary = json.load(open("config/parameters.json"))
    
    
    subset = get_subset(dictionary['read_1'],dictionary['read_2'],dictionary['subset_size'])
    fastqc(dictionary,subset)
    cutadapt(dictionary,subset)
    kallisto(dictionary,subset)
    gene_counts(dictionary,subset)
    
if __name__ == '__main__':
    main()