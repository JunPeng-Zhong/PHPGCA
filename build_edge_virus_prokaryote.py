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
###############################  virus-host   ##################################
################################################################################

# crispr
output_file = 'out/crispr_out.tab'
crispr_pred = {}
with open(output_file) as file_out:
    for line in file_out.readlines():
        parse = line.replace("\n", "").split("\t")
        virus = parse[0]
        prokaryote = parse[1].split('|')
        prokaryote = prokaryote[0].split('.')[0] if len(prokaryote) == 1 else prokaryote[1].split('.')[0]
        ident = float(parse[-3])
        length = float(parse[-2])
        slen = float(parse[-1])
        if virus not in crispr_pred and length/slen > 0.9 and ident > 0.9:
            crispr_pred[virus] = prokaryote

pkl.dump(crispr_pred, open('out/crispr_pred.dict', 'wb'))


# add connections between prokaryotes and viruses
tab_file_list = os.listdir(blast_tab_out)
prokaryote2virus = {}
for file in tab_file_list:
    virus_id_list = []
    with open(blast_tab_out+file) as file_in:
        for line in file_in.readlines():
            tmp = line.split('\t')
            virus_id = tmp[0].split('.')[0]
            prokaryote_id = tmp[1].split('.')[0]
            try:
                prokaryote2virus[prokaryote_id].append(virus_id)
            except:
                prokaryote2virus[prokaryote_id] = [virus_id]


# De-duplication
for key in prokaryote2virus:
    prokaryote2virus[key] = list(set(prokaryote2virus[key]))


# Save the virus-host graph
with open("out/phage_host.ntw", 'w') as file_out:
    for prokaryote in prokaryote2virus:
        for virus in prokaryote2virus[prokaryote]:
            _ = file_out.write(prokaryote + "," + virus + "\n")