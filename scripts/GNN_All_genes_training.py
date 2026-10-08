import pandas as pd
import numpy as np
from functions.CreateSTGraph import CreateSTGraph

from sklearn.preprocessing import RobustScaler
from torch_geometric.nn import GCNConv
import torch.nn as nn
from torchmetrics.classification import BinaryRecall,BinaryPrecision
from random import shuffle
import torch.optim as optim
import torch
from pathlib import Path
import pickle as pkl
from sklearn.metrics import precision_score, recall_score, f1_score, balanced_accuracy_score, roc_curve
import copy

metric_results = []
curve_results = []

target_genes = ["AMER1","BIRC6","CNTNAP2","CREBBP","CTNNA2","CTNND2","EP300","ERBB3","ERBB4","FAM47C","MALAT1","MTOR","MUC4","NF1","SMAD4","STAG2","TNC",'ACVR2A','AKAP9','APC','ARID1A','ATM','ATRX','BRAF','CASP8','CDH10','CDKN2A','CIC','CSMD3','CTNNB1','DCC','EGFR','FAM135B','FAT1','FAT3','FAT4','FBXW7','KEAP1','KMT2A','KMT2C','KMT2D','KRAS','LRP1B','MUC16','NBEA','NFE2L2','NOTCH1','NSD1','PDE4DIP','PIK3CA','PTEN','PTPRD','PTPRT','RB1','RNF213','RNF43','ROBO2','SETD2','SPOP','TP53','TRRAP','ZFHX3']

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(device)

## Create the GCN class which inherits from torch.nn.Module
class GCN(nn.Module):
  def __init__(self):
    super(GCN, self).__init__() ## calls constructor of the torch.nn.Module parent class
    self.conv1 = GCNConv(len(training_data.columns)-1, 1000)
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

def train(data):
    model.train()
    optimizer.zero_grad()
    out = model(data.x, data.edge_index)
    loss = criterion(out, data.y.view(-1, 1).float())
    loss.backward()
    optimizer.step()
    return loss

def evaluate(test_model, data):
    test_model.eval()
    with torch.no_grad():
        out = test_model(data.x, data.edge_index)
        loss = criterion(out, data.y.view(-1, 1).float())
        BinOut = torch.sigmoid(out)
        Recaller = BinaryRecall().to(device)
        Precisioner = BinaryPrecision().to(device)
        Recall = Recaller(BinOut,data.y.view(-1, 1).float())
        Precision = Precisioner(BinOut,data.y.view(-1, 1).float())
    return loss.item(), Precision, Recall

def select_cutoff(labels, scores):
    fpr, tpr, thresholds = roc_curve(labels, scores)

    valid = np.isfinite(thresholds)
    fpr = fpr[valid]
    tpr = tpr[valid]
    thresholds = thresholds[valid]

    best_index = np.argmax(tpr - fpr)

    return thresholds[best_index]


for target_gene in target_genes:
    model_path = 'data/large_model_RobustScale_scrambled.pt' ## Fill in path to the saved model (.pt file)
    print(target_gene)
    expr = pd.read_table('', header = 0, index_col=0, sep = ',') ## fill in path to the Combat batch corrected and unscaled expression training data for the target gene
    expr_test = pd.read_table('', header = 0, index_col=0, sep = ',') ## same as above but for the test data

    Y_train = expr['STATUS']
    Y_test = expr_test['STATUS']
    CancerType = expr['DISEASE']
    CancerType_test = expr_test['DISEASE']
    expr = expr.drop(columns=expr.columns[expr.columns.str.startswith("DISEASE")])
    expr = expr.drop(columns=['STATUS','n','w'])
    expr_test = expr_test.drop(columns=expr_test.columns[expr_test.columns.str.startswith("DISEASE")])
    expr_test = expr_test.drop(columns=['STATUS'])

    scaler = RobustScaler().fit(expr)
    training_data = pd.DataFrame(scaler.transform(expr), columns=expr.columns, index = expr.index)
    training_data = pd.concat([training_data,Y_train], axis=1)
    training_data
    testing_data = pd.DataFrame(scaler.transform(expr_test), columns=expr_test.columns, index = expr_test.index)
    testing_data = pd.concat([testing_data,Y_test], axis=1)
    with open(model_path+target_gene+'_RobustScaler.pkl','wb') as train_f:
        pkl.dump(scaler, train_f)

    ## Move the model and data to GPU
    model = GCN()
    model = model.to(device)

    criterion = nn.BCEWithLogitsLoss()
    optimizer = optim.Adam(model.parameters(), lr=0.00019971909638805173, weight_decay=0.01)

    def train(data):
        model.train()
        optimizer.zero_grad()
        out = model(data.x, data.edge_index)
        loss = criterion(out, data.y.view(-1, 1).float())
        loss.backward()
        optimizer.step()
        return loss

    def evaluate(test_model, data):
        test_model.eval()
        with torch.no_grad():
            out = test_model(data.x, data.edge_index)
            loss = criterion(out, data.y.view(-1, 1).float())
            BinOut = torch.sigmoid(out)
            Recaller = BinaryRecall().to(device)
            Precisioner = BinaryPrecision().to(device)
            Recall = Recaller(BinOut,data.y.view(-1, 1).float())
            Precision = Precisioner(BinOut,data.y.view(-1, 1).float())
        return loss.item(), Precision, Recall

    epochs = range(0, 100)
    train_losses = []
    test_losses = []
    test_precisions = []
    Recalls= []
    best_loss = float("inf")

    validation_graph = CreateSTGraph(
            testing_data,
            grid_size=32,
            num_circles=5,
            circle_radius=[4, 7],
            seed=1001)[0].to(device)

    for epoch in epochs:
        data = CreateSTGraph(training_data, grid_size=32, num_circles=5,circle_radius=[4,7], seed = epoch)[0]
        train_loss = train(data.to(device))
        data = data.to('cpu')
        train_losses.append(train_loss.to("cpu").detach().numpy())

        test_loss, test_precision, Recall = evaluate(model, validation_graph)
        test_losses.append(test_loss)
        test_precisions.append(test_precision.to('cpu'))
        Recalls.append(Recall.to('cpu').detach().numpy())
        # print(f"Epoch: {epoch}\tTrain Loss: {train_loss:.4f};\tTest_loss: {test_loss:.4f}", end = '\r')
        if epoch > 10 and test_loss < best_loss: 
            epoch_save = epoch
            best_loss = test_loss
            best_model_state = copy.deepcopy(model.state_dict())
            
            
    # Load the model with the lowest validation loss
    model.load_state_dict(best_model_state)

    model.eval()

    with torch.no_grad():
        logits = model(validation_graph.x,validation_graph.edge_index)
        test_scores = torch.sigmoid(logits).cpu().numpy().ravel()
        test_labels = validation_graph.y.cpu().numpy().astype(int).ravel()
        cutoff = select_cutoff(test_labels, test_scores)
        test_predictions = (test_scores >= cutoff).astype(int)

    precision = precision_score(
        test_labels,
        test_predictions,
        zero_division=0
    )

    recall = recall_score(
        test_labels,
        test_predictions,
        zero_division=0
    )

    f1 = f1_score(
        test_labels,
        test_predictions,
        zero_division=0
    )

    balanced_accuracy = balanced_accuracy_score(
        test_labels,
        test_predictions
    )

    # One summary row per gene
    metric_results.append({
        "gene": target_gene,
        "best_epoch": epoch_save,
        "validation_loss": best_loss,
        "cutoff": cutoff,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "balanced_accuracy": balanced_accuracy
    })

    # One row per test-graph spot
    curve_results.append(
        pd.DataFrame({
            "gene": target_gene,
            "spot": np.arange(len(test_labels)),
            "label": test_labels,
            "prediction_score": test_scores
        })
    )

    torch.save({"model_state_dict": model.state_dict(),
                "cutoff": cutoff},
                model_path + target_gene + "_combat_weight_nocnv.pt")
    print(
        f"Saved model from epoch {epoch_save}; "
        f"loss: {best_loss:.4f}; "
        f"precision: {precision:.4f}; "
        f"recall: {recall:.4f}; "
        f"F1: {f1:.4f}; "
        f"balanced accuracy: {balanced_accuracy:.4f}"
    )

output_path = '' ## Fill in path to the output directory where the CSV files will be saved

metric_results = pd.DataFrame(metric_results)

curve_results = pd.concat(
    curve_results,
    ignore_index=True
)

metric_results.to_csv(
    output_path + 'GCN_test_metrics_bad_GLM.csv',
    index=False
)

curve_results.to_csv(
    output_path + 'GCN_test_curve_data_bad_GLM.csv',
    index=False
)
# APC
# saved model from Epoch: 69       Train Loss: 0.1366      Test Loss: 0.1421      Test Precision: 0.9083  Test Recall: 0.9854
# ARID1A
# saved model from Epoch: 58       Train Loss: 0.2238      Test Loss: 0.1422      Test Precision: 0.9062  Test Recall: 0.9854
# ATRX
# saved model from Epoch: 56       Train Loss: 0.1857      Test Loss: 0.1671      Test Precision: 0.9259  Test Recall: 0.9709
# BRAF
# saved model from Epoch: 69       Train Loss: 0.3476      Test Loss: 0.3841      Test Precision: 0.9755  Test Recall: 0.9660
# CIC
# saved model from Epoch: 66       Train Loss: 0.2063      Test Loss: 0.1808      Test Precision: 0.9146  Test Recall: 0.9879
# CTNNB1
# saved model from Epoch: 52       Train Loss: 0.1236      Test Loss: 0.1309      Test Precision: 0.9161  Test Recall: 0.9806
# EGFR
# saved model from Epoch: 69       Train Loss: 0.4066      Test Loss: 0.3918      Test Precision: 0.9333  Test Recall: 0.9854
# KEAP1
# saved model from Epoch: 69       Train Loss: 0.3147      Test Loss: 0.3896      Test Precision: 0.9630  Test Recall: 0.9466
# KRAS
# saved model from Epoch: 69       Train Loss: 0.2746      Test Loss: 0.3272      Test Precision: 0.9485  Test Recall: 0.9393
# NFE2L2
# saved model from Epoch: 69       Train Loss: 0.2938      Test Loss: 0.2570      Test Precision: 0.9442  Test Recall: 0.9854
# NOTCH1
# saved model from Epoch: 54       Train Loss: 0.5825      Test Loss: 0.5630      Test Precision: 0.7432  Test Recall: 0.9976
# NSD1
# saved model from Epoch: 64       Train Loss: 0.1675      Test Loss: 0.1939      Test Precision: 0.9371  Test Recall: 0.9757
# PTEN
# saved model from Epoch: 69       Train Loss: 0.1075      Test Loss: 0.0992      Test Precision: 0.9441  Test Recall: 0.9830
# ROBO2
# saved model from Epoch: 69       Train Loss: 0.6089      Test Loss: 0.6057      Test Precision: 0.6338  Test Recall: 1.0000
# SETD2
# saved model from Epoch: 69       Train Loss: 0.3657      Test Loss: 0.3980      Test Precision: 0.9524  Test Recall: 0.9709
# TP53
# saved model from Epoch: 67       Train Loss: 0.1380      Test Loss: 0.1184      Test Precision: 0.9538  Test Recall: 0.9515
# ZFHX3
# saved model from Epoch: 66       Train Loss: 0.2002      Test Loss: 0.3314      Test Precision: 0.9131  Test Recall: 0.9442