import subprocess
import os
import argparse

parser = argparse.ArgumentParser(description='Switch the dataset')
parser.add_argument('--dataset', 
                    choices=['Clean', 'CHERRY', 'GPIC', 'HIC', 'case1', 'case2'], 
                    default='CHERRY', 
                    type=str, 
                    help='choose dataset')

inputs = parser.parse_args()

print('dataset: {0}'.format(inputs.dataset))

GCN_data = './GCN_data/'
dataset = './dataset/'
name_list = './name_list.csv'
prokaryote = './prokaryote/'
test_contigs = './test_contigs.fa'
out = './out/'
saved_model = './saved_model'
pred = './pred'
tmp_pred = './tmp_pred/'
edge_index = './edge_index/'
single_contig = './single_contig/'

change_list = [GCN_data, dataset, name_list, prokaryote, test_contigs, out, saved_model, pred, tmp_pred, edge_index, single_contig]

for fi in change_list:
    if os.path.exists(fi):
        print('{0} is exist and remove...'.format(fi))
        cmd = 'rm -rf {0}'.format(fi)
        subprocess.check_call(cmd, shell=True)

if(inputs.dataset == 'Clean'):
    for fi in ['Split_files', 'input', 'train_phage', 'all_proteins', 'allHOST.fasta']:
        if os.path.exists(fi):
            print('{0} is exist and remove...'.format(fi))
            cmd = 'rm -rf {0}'.format(fi)
            subprocess.check_call(cmd, shell=True)
    print('Only clean...')
    
else:
    print('Switch the dataset to {0}'.format(inputs.dataset))

    cmd = 'tar -zxvf datas/{0}.tar.gz'.format(inputs.dataset)
    subprocess.check_call(cmd, shell=True)

    os.makedirs(pred)
    os.makedirs(tmp_pred)

    print('Finish...')
