import os
import subprocess
import Bio
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

genome_list = os.listdir("prokaryote")

os.makedirs('CRISPR', exist_ok=True)
for genome in genome_list:
    CRISPR_cmd = 'java -cp ~/work/CRT1.2-CLI.jar crt prokaryote/'+ genome +' CRISPR/'+genome.split(".")[0]+'.cri'
    print("Capturing CRISPR")
    _ = subprocess.check_call(CRISPR_cmd, shell=True)


def special_match(strg, search=re.compile(r'[^ACGT]').search):
    return not bool(search(strg))

accession_file_list = os.listdir("CRISPR/")

os.makedirs('CRISPR_fasta')

for accession_file in accession_file_list:
    CRISPR_list = []
    with open('CRISPR/'+accession_file) as file_in:
        txt = list(file_in.readlines())
        for i in range(len(txt)):
            if 'No CRISPR elements were found.' in txt[i]:
                break
            try:
                CRISPR_list.append(txt[i].split('\t')[3])
            except:
                continue
    # remove nonCRISPR
    clean_CRISPR_list = []
    for seq in CRISPR_list:
        if special_match(seq) and seq != '':
            clean_CRISPR_list.append(seq)
    CRISPR_list = clean_CRISPR_list
    # save file
    if len(CRISPR_list) > 0:
        cnt = 1
        accession = accession_file.split('.')[0]
        record_list = []
        for CRISPR in CRISPR_list:
            rec = SeqRecord(Seq(CRISPR), id=accession+ "." + str(cnt), description="")
            record_list.append(rec)
            cnt += 1
        _ = SeqIO.write(record_list, "CRISPR_fasta/"+accession+".fasta", "fasta")

# 将不同的Accession的crispr整合到一个文件
crispr_list = os.listdir('CRISPR_fasta')
record_list = []
for acc in crispr_list:
    for record in SeqIO.parse('CRISPR_fasta/'+acc, 'fasta'):
        record_list.append(record)

SeqIO.write(record_list, 'allCRISPR.fasta','fasta') 

# 构造所有prokaryotes的crispr组成的database
make_blast_cmd = 'makeblastdb -in allCRISPR.fasta' + ' -dbtype nucl -parse_seqids -out dataset/crispr_db/allCRISPRs'
subprocess.check_call(make_blast_cmd, shell=True)