import os
import json
import glob
import numpy as np 
import pandas as pd


def get_subset(read1, read2, size):
    r1 = glob.glob(read1)
    r2 = glob.glob(read2)
    
    paired_end_reads = []
    for read in r2:
        sample = read.split("/")[-1].split("_")[0]
        if sample in "".join(r1):
            paired_end_reads.append(read)
        
    subset = paired_end_reads[:size]
    return subset

def fastqc(dictionary,subset):
    
    for sample in subset:
        s1 = sample.replace("_2.","_1.")
        s2 = sample
        command1 = "/opt/FastQC/fastqc " + s1 + " --outdir=" + dictionary['fast_out_path']
        command2 = "/opt/FastQC/fastqc " + s2 + " --outdir=" + dictionary['fast_out_path']

        os.system(command1)
        os.system(command2)
    
    os.chdir(dictionary['fast_out_path'])

    if not os.path.isdir('unzipped'):
        os.system('mkdir ' + 'unzipped')

    os.system('unzip \*.zip -d unzipped')
    os.chdir('..')
    return
    
def check_fast_qc(dictionary):
    #Parse and check for bad samples
    failed = []
    os.chdir('unzipped')

    for report in os.listdir():
        if report.startswith("SRR"):
            os.chdir(report)
            with open('summary.txt') as f:
                flag = f.readline()[:4]
                if flag != 'PASS':
                    failed.append(report)
            os.chdir('..')
        else:
            continue
    return failed

def cutadapt(dictionary,subset):
    os.chdir('cutadapt')

    for sample in subset:
        s1 = sample.replace("_2.","_1.")
        f1 = s1.split('/')[-1]

        s2 = sample
        f2 = s2.split('/')[-1]

        command1 = f"cutadapt -j 4 -a {dictionary['adapter_sequence']} -o {f1} ../../../{s1}"
        command2 = f"cutadapt -j 4 -a {dictionary['adapter_sequence']} -o {f2} ../../../{s2}"

        os.system(command1)
        os.system(command2)
        
    os.chdir('..')
    return

def kallisto(dictionary,subset):
    os.chdir('kallisto')
    
    for sample in subset:
        s1 = sample.replace("_2.","_1.")    
        s2 = sample

        outpath = s1.split('/')[-1].split('_')[0]

        if not os.path.isdir(outpath):
            os.system('mkdir ' + outpath)

        command = f"/opt/kallisto_linux-v0.42.4/kallisto quant -i {dictionary['idx']} -o {outpath} -b {dictionary['num_bootstraps']} ../../../{s1} ../../../{s2}"

        os.system(command)
    os.chdir('..')
    return
    
def gene_counts(dictionary,subset):
    
    tsvs = glob.glob("kallisto/SRR*/*.tsv")
    print(tsvs)
    gene_counts = pd.read_csv(tsvs[0],sep='\t')
    gene_counts = gene_counts[['target_id','est_counts']].rename(columns={'est_counts': tsvs[0].split('/')[-2]})
    
    for tsv in tsvs[1:]:
        df = pd.read_csv(tsv,sep='\t')
        df = df[['est_counts','target_id']].rename(columns={'est_counts': tsv.split('/')[-2]})
        gene_counts = gene_counts.join(df.set_index('target_id'), on='target_id')
        
    gene_counts.to_csv('kallisto/gene_counts.csv',index=False)
    return