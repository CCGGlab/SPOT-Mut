import matplotlib as plt
import pandas as pd
import os
import torch
import pickle as pkl

base_path = "data/attributions/"
top_k = 25
dataset = pd.read_csv(os.path.join("data", "LogNormalized_ALL_counts.csv.gz"), index_col=0).T
with open('data/TCGA_test_data_graphs_1.pkl','rb') as test_f:
    test_data_p1 = pkl.load(test_f)
with open('data/TCGA_test_data_graphs_2.pkl','rb') as test_f:
    test_data_p2 = pkl.load(test_f)
test_data = {**test_data_p1, **test_data_p2}
# ---- Get feature names ----
names = dataset.columns[:-1].tolist()

# ---- Collect tensors ----
def Get_Global_means(base_path, pkl_suffix):
    all_means = []

    for project in test_data.keys():
        pkl_path = os.path.join(base_path, f"{project}_rows_{pkl_suffix}.pkl")
        col_means = torch.load(pkl_path)   # shape: (13000,)
        all_means.append(col_means)

    # ---- Stack and average across projects ----
    all_means = torch.stack(all_means)          # shape: (num_projects, 13000)
    return(torch.mean(all_means, dim=0))

avg_attribution = Get_Global_means(base_path, 'abs')
# ---- Sort in descending order ----
sorted_vals, sorted_idx = torch.sort(avg_attribution, descending=True)

# ---- Extract top-k ----
top_vals = sorted_vals[:top_k].cpu().numpy()
top_idx = sorted_idx[:top_k].cpu().numpy()
top_names = [names[i] for i in top_idx]

# ---- define colors ----
avg_pos_attr = Get_Global_means(base_path, 'pos').to('cpu')
avg_neg_attr = Get_Global_means(base_path, 'neg').to('cpu')
geneCategories = {"Apoptosis":["BAX","FAS","AEN","TRIAP1","RPS27L","TNFRSF10B","TNFRSF10C","PHLDA3"],
                  "Metabolism":["FDXR","C1QBP"],
                  'DNA damage response': ["DDB2","POLH","PTCHD4"],
                  "p53 Autoregulation": ["CDKN2A","MDM2","RPL22L1"],
                  "Other": ["TM7SF3","PHLDA3","SPATA18","EDA2R","SHPK", "CDKN1A"]}

# ---- Category colors ----
category_colors = {
    "Apoptosis": "#D55E00",
    "Metabolism": "#009E73",
    "DNA damage response": "#0072B2",
    "p53 Autoregulation": "#CC79A7",
    "Other": "black",
    "Not Interacting": "#999999"
}

# ---- Create gene -> category mapping ----
gene_to_category = {}

for category, genes in geneCategories.items():
    for gene in genes:
        if gene not in gene_to_category:   # keeps first category if duplicated
            gene_to_category[gene] = category

# Category for top genes
top_categories = [
    gene_to_category.get(gene, "Not Interacting")
    for gene in top_names
]

# Label colors
label_colors = [
    category_colors[cat]
    for cat in top_categories
]

# ---- Plot ----
plt.figure(figsize=(5, 8))
plt.rcParams['font.family'] = 'Arial'

plt.barh(range(top_k), top_vals, color=label_colors)
ytick_addon = [" (-)" if avg_neg_attr[i] >= avg_pos_attr[i] else ' (+)' for i in top_idx]
top_names_addon = [top_names[i] + ytick_addon[i] for i in range(len(ytick_addon))]
plt.yticks(range(top_k), top_names_addon)

# Color tick labels
ax = plt.gca()

# for tick in ax.get_yticklabels():
#     tick.set_fontweight('bold')

ax.invert_yaxis()

plt.xlabel("Average Absolute Attribution")


for i, v in enumerate(top_vals):
    plt.text(v, i, f"{v:.5f}",
             va='center', ha='right',
             color='white')

import matplotlib.ticker as ticker
ax.xaxis.set_major_formatter(
    ticker.ScalarFormatter(useMathText=True)
)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

# ---- Legends ----
from matplotlib.patches import Patch

category_legend = [
    Patch(facecolor=c, label=cat)
    for cat, c in category_colors.items()
]

plt.legend(
    handles=category_legend,
    loc="lower right",
    frameon=False
)

plt.tight_layout()
# plt.savefig('guided_backprop_pathway.svg', transparent=True, facecolor='none', format = 'svg')
plt.show()