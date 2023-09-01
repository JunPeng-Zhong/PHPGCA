import os
import subprocess
import Bio
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

os.makedirs('./dataset/virus_db', exist_ok=True)

# 构造所有prokaryotes的crispr组成的database
make_blast_cmd = 'makeblastdb -in dataset/nucl.fasta' + ' -dbtype nucl -parse_seqids -out dataset/virus_db/allVIRUS'
subprocess.check_call(make_blast_cmd, shell=True)