import os
import json
import glob
import numpy as np 
import pandas as pd
import subprocess


def create_folders():
    if not os.path.isdir('data'):
        os.system('mkdir data')
    if not os.path.isdir('data/fastqc_results'):
        os.system('mkdir data/fastqc_results')
    if not os.path.isdir('data/cutadapt'):
        os.system('mkdir data/cutadapt')
    if not os.path.isdir('data/kallisto'):
        os.system('mkdir data/kallisto')
    return

def clean():
    if os.path.isdir('data'):
        os.system('rm -R data')
    return

def get_subset(read1, read2, size):
    r1 = glob.glob(read1)
    r2 = glob.glob(read2)
    
    paired_end_reads = []
    for read in r2:
        sample = read.split("/")[-1].split("_")[0]
        if sample in "".join(r1):
            paired_end_reads.append(read)
        
    subset = paired_end_reads[100:117]
    return subset

def fastqc(dictionary,subset):
    
    for sample in subset:
        s1 = sample.replace("_2.","_1.")
        s2 = sample
        command1 = "/opt/FastQC/fastqc " + s1 + " --outdir=" + dictionary['fast_out_path']
        command2 = "/opt/FastQC/fastqc " + s2 + " --outdir=" + dictionary['fast_out_path']

        os.system(command1)
        os.system(command2)

    zip_files = glob.glob(dictionary['fast_out_path']+"/*.zip")
    
    for file in zip_files:
        os.system(f"unzip {file} -d {dictionary['fast_out_path']}/unzipped")
        
    return
    
def check_fast_qc(dictionary):
    #Parse and check for bad samples
    failed = []
    
    for report in os.listdir('unzipped'):
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

    for sample in subset:
        s1 = sample.replace("_2.","_1.")
        f1 = s1.split('/')[-1]

        s2 = sample
        f2 = s2.split('/')[-1]

        command1 = f"cutadapt -j 4 -a {dictionary['adapter_sequence']} -o data/cutadapt/{f1} {s1}"
        command2 = f"cutadapt -j 4 -a {dictionary['adapter_sequence']} -o data/cutadapt/{f2} {s2}"

        os.system(command1)
        os.system(command2)
        
    return

def kallisto(dictionary,subset):
    
    for sample in subset:
        s1 = sample.replace("_2.","_1.")    
        s2 = sample

        outpath = "data/kallisto/" + s1.split('/')[-1].split('_')[0]

        if not os.path.isdir(outpath):
            os.system('mkdir ' + outpath)

        command = f"/opt/kallisto_linux-v0.42.4/kallisto quant -i {dictionary['idx']} -o {outpath} -b {dictionary['num_bootstraps']} -t {dictionary['num_bootstraps']} {s1} {s2}"

        os.system(command)
    return
    
def gene_counts(dictionary,subset):
    
    tsvs = glob.glob("data/kallisto/SRR*/*.tsv")
    gene_counts = pd.read_csv(tsvs[0],sep='\t')
    gene_counts = gene_counts[['target_id','est_counts']].rename(columns={'est_counts': tsvs[0].split('/')[-2]})
    
    for tsv in tsvs[1:]:
        df = pd.read_csv(tsv,sep='\t')
        df = df[['est_counts','target_id']].rename(columns={'est_counts': tsv.split('/')[-2]})
        gene_counts = gene_counts.join(df.set_index('target_id'), on='target_id')
            
    names = pd.read_csv('/datasets/srp073813/reference/Gene_Naming.csv')
    names = names.loc[(names['chr'] != 'chrY') & (names['chr'] != 'chrX')]
    
    gene_counts = gene_counts.loc[gene_counts['target_id'].str.startswith('NM_')]
    joined = pd.merge(names, gene_counts, how='left', left_on='refseq', right_on='target_id')
    
    filtered = joined.filter(regex='^(SRR*|refseq$)', axis=1).dropna()
    filtered.loc[:,filtered.columns[1:]] = filtered.filter(regex='^(SRR*)', axis=1).round().astype('int64') + 1
    filtered.to_csv('data/filtered_gene_counts.csv',index=False)
    
    return

def pca():
    subprocess.call(['/opt/conda/bin/Rscript',  '--vanilla', 'pca.r'])
    return 