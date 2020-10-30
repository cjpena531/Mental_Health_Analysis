import os 
import json

lst = os.listdir('/datasets/srp073813')

def get_names(lst):
    with open('names.txt', 'w') as f:
        for item in lst:
            f.write('%s\n' % item)
