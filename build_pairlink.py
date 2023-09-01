import time
import os
import sys
import Bio
import logging
import argparse
import subprocess
import scipy as sp
import numpy as np
import pandas as pd
import pickle as pkl
import networkx as nx
import scipy.stats as stats
import scipy.sparse as sparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline


parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--mode', type=str, default = 'virus')
parser.add_argument('--threads', type=int, default=4, help='num_threads for blastn')
parser.add_argument('--refresh', type=str, default='False', help='Refresh the blastn')
parser.add_argument('--refresh-allHOST', type=str, default='False', help='Refresh the allHOST.fasta')

inputs = parser.parse_args()


# Defined folder
prokaryote_fn= "prokaryote/"
new_prokaryote = "new_prokaryote/"
blast_database_out = "blast_db1/"
new_blast_database_out = "new_blast_db/"
blast_tab_out = "blast_tab1/"
new_blast_tab_out = "new_blast_tab/"
Knowledge_graph = "Cyber_data/"

################################################################################
############################  Check the folder #################################
################################################################################
def check_folder(file_name):
    if not os.path.exists(file_name):
        _ = os.makedirs(file_name)
    else:
        print("folder {0} exist... cleaning dictionary".format(file_name))
        if os.listdir(file_name):
            try:
                _ = subprocess.check_call("rm -rf {0}".format(file_name), shell=True)
                _ = os.makedirs(file_name)
                print("Dictionary cleaned")
            except:
                print("Cannot clean your folder... permission denied")
                exit(1)

if(inputs.refresh == 'True'):
    check_folder(blast_database_out)
    check_folder(blast_tab_out)
    if inputs.mode != 'virus':
        check_folder(new_blast_tab_out)
        check_folder(new_blast_database_out)



# combine phage file 
if(not os.path.exists('out/query.fa')):
    _ = subprocess.check_call("cat dataset/nucl.fasta single_contig/* > out/query.fa", shell=True)




################################################################################
###############################  Run CRISPR   ##################################
################################################################################

print('Runing CRISPR')
query_file = "out/test.fa"
db_host_crispr_prefix = "dataset/crispr_db/allCRISPRs"
output_file = "out/crispr_out.tab"
crispr_call = NcbiblastnCommandline(query=query_file,db=db_host_crispr_prefix,out=output_file,outfmt="6 qseqid sseqid evalue pident length slen", evalue=1,gapopen=10,penalty=-1,
                                  gapextend=2,word_size=7,dust='no',
                                 task='blastn-short',perc_identity=90,num_threads=16)
crispr_call()


# crispr_pred = {}
# with open(output_file) as file_out:
#     for line in file_out.readlines():
#         parse = line.replace("\n", "").split("\t")
#         virus = parse[0]
#         prokaryote = parse[1]
#         prokaryote = prokaryote.split('.')[0]
#         ident = float(parse[-3])
#         length = float(parse[-2])
#         slen = float(parse[-1])
#         if virus not in crispr_pred and length/slen > 0.9 and ident > 0.9:
#             crispr_pred[virus] = prokaryote
# 
# pkl.dump(crispr_pred, open('out/crispr_pred.dict', 'wb'))



################################################################################
###############################  Speed UP ######################################
################################################################################
print('Building allHOST.fasta')

if(inputs.refresh_allHOST == 'True' and os.path.exists('allHOST.fasta')):
	cmd = 'rm allHOST.fasta'
	print('Deleting allHOST.fasta')
	subprocess.check_call(cmd, shell=True)


if(not os.path.exists('allHOST.fasta')):
    proka_dir = 'prokaryote/'
    proka_list = os.listdir(proka_dir)

    records = []
    for proka in proka_list:
        cnt = 1
        proka_name = proka.split('.')[0]
        for record in SeqIO.parse(proka_dir+proka, 'fasta'):
            record.id = proka_name+'.'+str(cnt)
            cnt += 1
            records.append(record)
    SeqIO.write(records, 'allHOST.fasta', 'fasta')

else:
    print('allHOST.fasta exist...')

################################################################################
###############################  Run BLASTN(train)   ###########################
################################################################################

print('Runing BLASTN')
make_blast_cmd = 'makeblastdb -in allHOST.fasta -dbtype nucl -parse_seqids -out '+ blast_database_out +'allHOST'
print("Creating blast database...")
_ = subprocess.check_call(make_blast_cmd, shell=True)
blast_cmd = 'blastn -query out/query.fa -db blast_db1/allHOST -outfmt 6 -out '+ blast_tab_out + 'allHOST.tab -num_threads {0} -mt_mode 1'.format(inputs.threads)
print("Running blastn...")
start = time.time()
_ = subprocess.check_call(blast_cmd, shell=True)
end = time.time()
print('Blastn use {0:.2f}s'.format(end-start))
    


################################################################################
###############################  Run BLASTN(test)   ############################
################################################################################

if inputs.mode != 'virus':
    genome_list = os.listdir(new_prokaryote)
    for genome in genome_list:
        make_blast_cmd = 'makeblastdb -in '+ new_prokaryote + genome +' -dbtype nucl -parse_seqids -out '+ new_blast_database_out + genome.split(".")[0]
        print("Creating blast database...")
        _ = subprocess.check_call(make_blast_cmd, shell=True)
        blast_cmd = 'blastn -query out/query.fa -db new_blast_db/'+genome.split(".")[0]+' -outfmt 6 -out '+ new_blast_tab_out + genome.split(".")[0]+'.tab -num_threads 32 -mt_mode 1'
        print("Running blastn...")
        _ = subprocess.check_call(blast_cmd, shell=True)
