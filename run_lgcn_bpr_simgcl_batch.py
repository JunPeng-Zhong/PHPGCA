import numpy as np
import pandas as pd
import torch
import pickle as pkl
import os
from torch_geometric.utils import from_scipy_sparse_matrix
from torch_sparse import SparseTensor
from sklearn.preprocessing import MinMaxScaler
import utils
import random
from torch_geometric.nn import GraphSAGE
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import RobustScaler, StandardScaler
import argparse
import model2
import torch.nn.functional as F
from torch.utils.data import DataLoader
from Bio import SeqIO
import time

parser = argparse.ArgumentParser(description='run pyg')
parser.add_argument('--device', default='cpu', type=str, help='Device of Training')
parser.add_argument('--epochs', default=4000, type=int, help='Epochs to train')
parser.add_argument('--train-mode', default='retrain', type=str, help='Mode for training')
parser.add_argument('--threshold', default=0.9, type=float, help='ACC threshold of training')
parser.add_argument('--lr', default=0.001, type=float, help='Learning rate of training')
parser.add_argument('--mode', default='virus', type=str, help='Mode of Training')
parser.add_argument('--weight-decay', default=0, type=float, help='Weight Decay of optimizer')
parser.add_argument('--cl-rate', default=0.5, type=float, help='CL rate for CL loss')
parser.add_argument('--eps', default=0.1, type=float, help='eps for CL')
parser.add_argument('--hidden-channels', default=256, type=int, help='Hidden Channels For training')
parser.add_argument('--num-layers', default=3, type=int, help='Layers of LightGCN')
parser.add_argument('--batch-size', default=512, type=int, help='Batch size of Training')
parser.add_argument('--topk', default=1, type=int, help='topk')

inputs = parser.parse_args()

seed = 42
np.random.seed(seed)
torch.random.manual_seed(seed)

device = torch.device(inputs.device)

# load data
adj         = pkl.load(open("GCN_data/graph.list",'rb'))
id2node     = pkl.load(open("GCN_data/id2node.dict",'rb'))
node2id     = pkl.load(open("GCN_data/node2id.dict", "rb" ))
idx_test    = pkl.load(open("GCN_data/test_id.dict", 'rb'))
node2label  = pkl.load(open("GCN_data/node2label.dict",'rb'))
crispr_pred = pkl.load(open('out/crispr_pred.dict', 'rb'))
prokaryote_df = pd.read_csv('dataset/prokaryote.csv')

# num test_virus
num_test = 0
for _ in SeqIO.parse('test_contigs.fa', 'fasta'):
    num_test += 1


# prokaryotes in the training set
trainable_host = []
for file in os.listdir('prokaryote/'):
    trainable_host.append(file.rsplit('.', 1)[0])


host2id = {}
label2hostid =  {}
trainable_host_idx = []
trainable_label = []
for idx, node in id2node.items():
    # if prokaryote
    if node in trainable_host:
        host2id[node] = idx
        trainable_host_idx.append(idx)
        trainable_label.append(node2label[node])
        if(node2label[node] not in label2hostid):
            label2hostid[node2label[node]] = []
        label2hostid[node2label[node]].append(idx)

trainable_label = list(set(trainable_label))

# Test preprocessing
virus2spe = {}
df_virus = pd.read_csv('dataset/virus.csv')

for ii in range(len(df_virus)):
    virus = df_virus.loc[ii]['Accession']
    spe = df_virus.loc[ii]['Species']

    virus2spe[virus] = spe

cherry2name = {}
name2cherry = {}
df_name = pd.read_csv('name_list.csv')

for ii in range(len(df_name)):
    name = df_name.loc[ii]['contig_name'].split('.')[0]
    cherry = df_name.loc[ii]['idx']

    cherry2name[cherry] = name
    name2cherry[name] = cherry

# GCN Model Data
edge_index = from_scipy_sparse_matrix(adj)[0]
edge_index = edge_index
adj_t = SparseTensor(row=edge_index[0], col=edge_index[1]).t().to(device)
# x = torch.tensor(features, dtype=torch.float, device=device)
x = None

# Pos & Neg preprocessing
preprocess_dir = 'edge_index/'
if(os.path.exists(preprocess_dir)):
    print(preprocess_dir+' exist and load data...')
    edge_index_pos_list = torch.load(open(preprocess_dir+'edge_index_pos_list', 'rb'))
    edge_index_pos = torch.load(open(preprocess_dir+'edge_index_pos', 'rb'))
    trainable_host_idx_set = torch.load(open(preprocess_dir+'trainable_host_idx_set', 'rb'))
    virus2pro_neg = torch.load(open(preprocess_dir+'virus2pro_neg', 'rb'))

else:
    os.makedirs(preprocess_dir)
    print('Preprocessing the edge_index data')
    edge_index_pos_list = [[],[]]

    for label in trainable_label:
        host_idxs = label2hostid[label]
        for host_idx in host_idxs:
            for idx in range(len(node2id)):
                if(idx not in trainable_host_idx and node2label[id2node[idx]] == label):
                    edge_index_pos_list[0].append(idx)
                    edge_index_pos_list[1].append(host_idx)
    edge_index_pos = torch.tensor(edge_index_pos_list, dtype=torch.long)

    trainable_host_idx_set = set(trainable_host_idx)
    virus2pro_neg = {}  # virus: list(neg prokaryotes)

    for virus, pro in zip(edge_index_pos_list[0], edge_index_pos_list[1]):
        virus2pro_neg[virus] = list(trainable_host_idx_set - set([pro])) 

    torch.save(edge_index_pos_list, preprocess_dir+'edge_index_pos_list')
    torch.save(edge_index_pos, preprocess_dir+'edge_index_pos')
    torch.save(trainable_host_idx_set, preprocess_dir+'trainable_host_idx_set')
    torch.save(virus2pro_neg, preprocess_dir+'virus2pro_neg')

print('Finished preprocess..')

def neg_sample(edge_index_pos_batch):
    edge_index_neg = [[],[]]
    for virus in edge_index_pos_batch[0].cpu().tolist():
        fake_pro = random.choice(virus2pro_neg[virus])
        edge_index_neg[0].append(virus)
        edge_index_neg[1].append(fake_pro)
    
    return torch.tensor(edge_index_neg, dtype=torch.long, device=device)


def train_accuracy():
    with torch.no_grad():
        total = 0
        correct = 0
        for i in range(len(encode)):
            if idx_test[id2node[i]] != 0 or i in trainable_host_idx:
                continue
            virus_feature = encode[i]
            max_pred = 0
            pred_label = ""
            for label in trainable_label:
                prokaryote_feature = encode[label2hostid[label]]
                preds = decoder(virus_feature, prokaryote_feature)
                for pred in preds:
                    if pred > max_pred:
                        max_pred = pred
                        pred_label = label
            if pred_label == node2label[id2node[i]]:
                correct+=1
            total += 1
    return correct/total

def train_accuracy2():
    with torch.no_grad():
        total = 0
        correct = 0
        for i in range(len(encode)):
            # if idx_test[id2node[i]] != 0 or i in trainable_host_idx:
            if 'cherry' in id2node[i] or id2node[i] not in idx_test:
                continue
            virus_feature = encode[i]
            max_pred = 0
            pred_label = ""

            prokaryote_feature = encode[trainable_host_idx]
            preds = decoder(virus_feature, prokaryote_feature)
            pred_label = node2label[id2node[trainable_host_idx[preds.argmax().cpu().item()]]]

            if pred_label == node2label[id2node[i]]:
                correct+=1
            total += 1
    return correct/total


# CL Loss
def cal_cl_loss(edge_index_pos_batch):
    emb_view1 = net(x, adj_t, perturbed=True)
    emb_view2 = net(x, adj_t, perturbed=True)

    v_idx = edge_index_pos_batch[0]
    p_idx = edge_index_pos_batch[1]

    virus_cl_loss = utils.InfoNCE(emb_view1[v_idx], emb_view2[v_idx], 0.2)
    proka_cl_loss = utils.InfoNCE(emb_view1[p_idx], emb_view2[p_idx], 0.2)

    return virus_cl_loss + proka_cl_loss



net = model2.SimGCL_Encoder(len(node2id), inputs.hidden_channels, eps=inputs.eps, device=device, num_layers=inputs.num_layers).to(device)
decoder = model2.DotDecoder().to(device)

params = list(net.parameters()) + list(decoder.parameters())
optimizer = torch.optim.Adam(params, lr=inputs.lr, weight_decay=inputs.weight_decay)
loss_func = torch.nn.BCEWithLogitsLoss()
scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=100, gamma=0.8)

if(inputs.train_mode == 'pretrain'):
    net_dict = torch.load(f"saved_model/encoder.pkl", map_location='cpu')
    net.load_state_dict(net_dict)

def collate_fn(data):
    return torch.hstack(data)

batch_size = inputs.batch_size

edge_index_pos_datset = utils.PosEdge(edge_index_pos)
edge_index_pos_dataloader = DataLoader(edge_index_pos_datset, batch_size=batch_size, shuffle=True, collate_fn=collate_fn)


#################################################################
##########################  Training  ###########################
#################################################################

early_stop = 5
cur_test_acc = 0.0
EPS = 1e-15

if inputs.train_mode == 'retrain':
    _ = net.train()
    _ = decoder.train()
    for epoch in range(1, 1+inputs.epochs):

        # Batch
        batch_cal = utils.AverageMeter()
        for batch, edge_index_pos_batch in enumerate(edge_index_pos_dataloader):

            #print(epoch)
            encode = net(x, adj_t)
            
            edge_index_neg_batch = neg_sample(edge_index_pos_batch)

            loss = 0

            virus_feat = encode[edge_index_pos_batch[0]]
            pro_feat_pos = encode[edge_index_pos_batch[1]]
            pro_feat_neg = encode[edge_index_neg_batch[1]]

            # bpr loss
            pred_pos = decoder(virus_feat, pro_feat_pos)
            pred_neg = decoder(virus_feat, pro_feat_neg)
            bpr_loss = torch.mean(F.softplus(pred_neg - pred_pos)) # mean or sum

            cl_loss = inputs.cl_rate * cal_cl_loss(edge_index_pos_batch)

            loss = bpr_loss + cl_loss

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            batch_cal.update(loss.cpu().item())
        
        #print(loss.cpu().detach().numpy())
        if (epoch % 10 == 0):
            train_acc = train_accuracy2()
            print("[Epoch {0:>4d}] Train acc: {1:.2f}, Total Loss: {2:.3f}, BPR Loss: {3:.3f}, CL Loss: {4:.3f}".format(epoch, train_acc, loss.cpu().item(), bpr_loss.cpu().item(), cl_loss.cpu().item()))
            
            if train_acc > 0.8:
                break
                    

    torch.save(net.state_dict(), f'saved_model/encoder.pkl')
    torch.save(decoder.state_dict(), f'saved_model/decoder.pkl')



#################################################################
#########################  Prediction  ##########################
#################################################################

total_confident = 0
top_eq_conf = 0

# predicting host
if inputs.mode == 'virus':
    node2pred = {}
    with torch.no_grad():
        encode = net(x, adj_t)
        for i in range(len(encode)):
            confident_label = 'unknown'
            if id2node[i] not in idx_test or idx_test[id2node[i]] == 0:
                continue
            if id2node[i] in idx_test and idx_test[id2node[i]] == 1:
                confident_label = node2label[id2node[i]]
            virus_feature = encode[i]
            pred_label_score = []
            for label in set(trainable_label):
                if label == confident_label:
                    pred_label_score.append((label, 1))
                    continue
                prokaryote_feature = encode[label2hostid[label]]
                preds = decoder(virus_feature, prokaryote_feature)
                for pred in preds:
                    pred_label_score.append((label, torch.sigmoid(pred).detach().cpu().numpy()))
            node2pred[id2node[i]] = sorted(pred_label_score, key=lambda tup: tup[1], reverse=True)

            if(idx_test[id2node[i]] == 1 and confident_label in set(trainable_label)):
                total_confident += 1
                if(node2pred[id2node[i]][0][0] == confident_label):
                    top_eq_conf += 1

        for virus in crispr_pred:
            if virus not in node2pred:
                pred = prokaryote_df[prokaryote_df['Accession'] == crispr_pred[virus]]['Species'].values[0]
                node2pred[virus] = [(pred, 1)]

k = inputs.topk
data = {
    'contig': []
}

for i in range(k):
    data['top_{0}'.format(i+1)] = []

for node, pred in node2pred.items():
    data['contig'].append(cherry2name[node])
    for i in range(k):
        data['top_{0}'.format(i+1)].append(pred[i][0])

df_pred = pd.DataFrame(data=data)
df_pred.to_csv('final_prediction.csv', index=False)

print('Prediction Finished...')