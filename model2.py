import  torch
from    torch.nn import functional as F
# PYG GCN
from torch_geometric.nn.conv import LGConv


class DotDecoder(torch.nn.Module):
    def __init__(self):
        super().__init__()

    def forward(self, x1, x2):
        out = (x1*x2).sum(dim=-1)
        return out.reshape(-1)

#########################################################################
######################## SimGCL #########################################
#########################################################################
class SimGCL_Encoder(torch.nn.Module):
    def __init__(self, num_node, num_feat, eps, device, num_layers=3):
        super().__init__()
        self.eps = eps
        self.num_layers = num_layers
        self.device = device
        self.layers = torch.nn.ModuleList()

        # Init Embedding
        self.emb = torch.nn.Embedding(num_node, num_feat)
        torch.nn.init.xavier_normal_(self.emb.weight)

        for ii in range(num_layers):
            self.layers.append(LGConv())

    def forward(self, x2, adj_t, perturbed=False):
        all_embeddings = []

        # x1 = self.emb.weight
        x_out = self.emb.weight
        for layer in self.layers:
            x_out = layer(x_out, adj_t)
            if(perturbed):
                random_noise = torch.rand_like(x_out, device=self.device)
                x_out = x_out + torch.sign(x_out) * F.normalize(random_noise, dim=-1) * self.eps
            all_embeddings.append(x_out)
        
        all_embeddings = torch.stack(all_embeddings, dim=1)
        all_embeddings = torch.mean(all_embeddings, dim=1)

        return all_embeddings
    