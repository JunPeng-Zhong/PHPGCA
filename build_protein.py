import os
import sys
import Bio
import subprocess
import numpy as np
import pandas as pd
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
from sklearn.model_selection import KFold

######################################
######## INPUT PARAMETERS ############
######################################
parser = argparse.ArgumentParser(description='Build Protein for training virus')
parser.add_argument('--addition-test', default='False', help='Have additional test virus')

inputs = parser.parse_args()


custom_dir = 'custom_dataset/'
virus_dir = 'custom_dataset/all_virus/'
dataset_dir = 'dataset/'
os.makedirs(dataset_dir, exist_ok=True)



if(inputs.addition_test == 'False'):
    # Split All virus to train and test
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    test_idx = list(kf.split(range(len(os.listdir(virus_dir)))))[0][1]


    record_list_tra = []
    record_list_test = []

    virus_names = os.listdir(virus_dir)
    for idx, virus_file in enumerate(virus_names):
        record = SeqIO.read(virus_dir+virus_file, 'fasta')
        if(idx in test_idx):
            record_list_test.append(record)
        else:
            record_list_tra.append(record)


    SeqIO.write(record_list_tra, dataset_dir+'nucl.fasta', 'fasta')
    SeqIO.write(record_list_test, 'test_contigs.fa', 'fasta')


# Refresh, get Protein of Training virus 
prodigal_cmd = 'prodigal -i {0}nucl.fasta -a {0}protein.fasta -f gff -p meta'.format(dataset_dir)
print("Running prodigal...")
_ = subprocess.check_call(prodigal_cmd, shell=True)

proteins = []
contigs  = []
keywords = []
for record in SeqIO.parse('{0}protein.fasta'.format(dataset_dir), 'fasta'):
    name = record.id
    contigs.append(name.rsplit("_", 1)[0])
    proteins.append(name)
    keywords.append('hypothetical protein')

gene2genome_df = pd.DataFrame({'protein_id': proteins, 'contig_id': contigs, 'keywords': keywords})
gene2genome_df.to_csv('{0}database_gene_to_genome.csv'.format(dataset_dir), index=False)