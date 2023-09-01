from Bio import SeqIO
import subprocess
import utils

#####################################################################
#########################   processing    ###########################
#####################################################################

utils.check_folder("input")
utils.check_folder("pred")
utils.check_folder('Split_files')
utils.check_folder('tmp_pred')
utils.check_folder("train_phage")

for record in SeqIO.parse('dataset/nucl.fasta', 'fasta'):
    _ = SeqIO.write(record, 'train_phage/'+record.id, 'fasta')

#####################################################################
#########################  Start Program  ###########################
#####################################################################
# split into sub files
cnt = 0
file_id = 0
records = []
for record in SeqIO.parse('test_contigs.fa', 'fasta'):
    if cnt !=0 and cnt%2000 == 0:
        SeqIO.write(records, f"Split_files/contig_{file_id}.fasta","fasta") 
        records = []
        file_id+=1
        cnt = 0
    seq = str(record.seq)
    seq = seq.upper()
    if len(record.seq) > 0:
        records.append(record)
        cnt+=1

SeqIO.write(records, f"Split_files/contig_{file_id}.fasta","fasta")
file_id+=1


# run sub files
for i in range(file_id):
    cmd = f"mv Split_files/contig_{i}.fasta input/"
    try:
        out = subprocess.check_call(cmd, shell=True)
    except:
        print(f"Moving file Error for file contig_{i}")
        exit()