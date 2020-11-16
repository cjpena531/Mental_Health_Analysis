import os
import json
import glob


def get_subset(read1, read2, size):
    r1 = glob.glob(read1)
    r2 = glob.glob(read2)
    
    paired_end_reads = []
    for read in r2:
        sample = read.split("/")[3].split("_")[0]
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

        command1 = f"cutadapt -j 4 -a {dictionary['adapter_sequence']} -o {f1} {s1}"
        command2 = f"cutadapt -j 4 -a {dictionary['adapter_sequence']} -o {f2} {s2}"

        os.system(command1)
        os.system(command2)
        
    os.chdir('..')

def kallisto(dictionary,subset):
    os.chdir('kallisto')
    
    for sample in subset:
        s1 = sample.replace("_2.","_1.")    
        s2 = sample

        outpath = s1.split('/')[-1].split('_')[0]

        if not os.path.isdir(outpath):
            os.system('mkdir ' + outpath)

        command = f"/opt/kallisto_linux-v0.42.4/kallisto quant -i {dictionary['idx']} -o {outpath} -b {dictionary['num_bootstraps']} {s1} {s2}"

        os.system(command)
    os.chdir('..')