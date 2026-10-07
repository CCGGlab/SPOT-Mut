## This is the GCN model. Use the weights from the saved model to make predictions on new data.
## Load the for any gene other than TP53 with 
## model = GCN(n_features = len(gene_set))
## saved_model = torch.load(gene+'_combat_weight_nocnv.pt', weights_only=False)
## model.load_state_dict(saved_model["model_state_dict"])
## OR
## If you are loading the TP53 model from the data/TP53_model.pt file, use the following code:
## model = torch.load('data/TP53_model.pt', weights_only=False)

## models can be downloaded from https://ccgg.ugent.be/shiny/spot-mut/, just select the gene of interest and click on the 'download model weights and scaler' button. The model weights and data scaler will be downloaded.

import torch.nn as nn
from torch_geometric.nn import GCNConv
class GCN(nn.Module):
  def __init__(self):
    super(GCN, self).__init__() ## calls constructor of the torch.nn.Module parent class
    self.conv1 = GCNConv(13013, 1000)
    self.conv3 = GCNConv(1000, 500)
    self.linear1 = nn.Linear(500, 256)
    self.linear2 = nn.Linear(256, 128)
    self.linear3 = nn.Linear(128, 64)
    self.linear4 = nn.Linear(64, 16)
    self.classifier = nn.Linear(16,1)  ## binary classification
    self.drop1 = nn.Dropout(0.5)
    self.drop3 = nn.Dropout(0.20)
    self.drop4 = nn.Dropout(0.5)
  def forward(self, x, edge_index):
    h = self.conv1(x, edge_index)
    h = h.relu()
    h = self.drop1(h)
    h = self.conv3(h, edge_index)
    h = h.relu()
    h = self.linear1(h)
    h = h.relu()
    h = self.drop3(h)
    h = self.linear2(h)
    h = h.relu()
    h = self.drop4(h)
    h = self.linear3(h)
    h = h.relu()
    h = self.linear4(h)
    h = h.relu()
    out = self.classifier(h)
    return out