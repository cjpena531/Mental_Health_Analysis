import os 
import json
import sys
sys.path.insert(0,sys.path[0] +'/src')
from etl import *

def main():
       
    dictionary = json.load(open("config/default_parameters.json"))
    subset = get_subset(dictionary['read_1'],dictionary['read_2'],dictionary['subset_size'])
    fastqc(dictionary,subset)
    cutadapt(dictionary,subset)
    kallisto(dictionary,subset)
    
if __name__ == '__main__':
    main()