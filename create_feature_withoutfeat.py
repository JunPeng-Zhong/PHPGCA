import os
import numpy as np
import Bio
from Bio import SeqIO
import pickle as pkl
import argparse



parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--mode', type=str, default = 'virus')
inputs = parser.parse_args()

os.makedirs('node_feature/', exist_ok=True)

def node2id(file_in_fn):
    # search files
    file_list = os.listdir(file_in_fn)
    file2idx = {}
    # convert to words
    for idx, file in enumerate(file_list):
        file2idx[file.rsplit('.', 1)[0]] = idx
    return file2idx


virus2id = node2id('train_phage/')
pkl.dump(virus2id, open('node_feature/virus.dict', 'wb'))

prokaryote2id = node2id('prokaryote/')
pkl.dump(prokaryote2id, open('node_feature/prokaryote.dict', 'wb'))

if inputs.mode == 'virus':
    test_virus2id = node2id('single_contig/')
    pkl.dump(test_virus2id, open('node_feature/test_virus.dict', 'wb'))
elif inputs.mode == 'prokaryote':
    test_prokaryote2id = node2id('new_prokaryote/')
    pkl.dump(test_prokaryote2id, open('node_feature/test_prokaryote.dict', 'wb'))
    test_virus2id = node2id('single_contig/')
    pkl.dump(test_virus2id, open('node_feature/test_virus.dict', 'wb'))
else:
    print('wrong parameters')
    exit()



