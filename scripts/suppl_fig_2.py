import matplotlib as plt
import pickle as pkl
import numpy as np

####panel a ####
with open('data/AUC_Final_C2_KEGG.pkl','rb') as acc_f:
    Accuracies = pkl.load(acc_f)
with open('data/r_auc_Final_C2_KEGG.pkl','rb') as r_acc_f:
    r_accuracies = pkl.load(r_acc_f)
# Sorting the bars based on the second item in the list (sorting_value)
means = {key: np.mean(values) for key, values in Accuracies.items()}
r_means = {key: np.mean(values) for key, values in r_accuracies.items()}
means = {key:means[key] - r_means[key] for key in means.keys()}
## filtering
# Top 10 positive
top_pos = sorted(
    [(k, v) for k, v in means.items() if v > 0],
    key=lambda x: x[1],
    reverse=True
)[:20]

# Top 10 negative
top_neg = sorted(
    [(k, v) for k, v in means.items() if v < 0],
    key=lambda x: x[1]
)[:20]

# Combine and sort for plotting
selected = dict(top_pos + top_neg)

sorted_acc = dict(
    sorted(selected.items(), key=lambda item: item[1], reverse=True)
)

# Extracting sorted x and y values
sorted_x_values = [key for key in sorted_acc.keys()]
sorted_y_values = [means[key] for key in sorted_acc.keys()]
bar_colors = ['#9E0142' if i <0 else '#5E4FA2' for i in sorted_y_values]
ytickz = [' '.join(hallmark.split('_')[1:]) for hallmark in [Hallmark[:-4] for Hallmark in sorted_x_values]]
Fontsize = 10
plt.rcParams['font.family'] = 'Arial'
plt.figure(figsize=(5, 10))
plt.barh(sorted_x_values, sorted_y_values, 0.4, color = bar_colors)
plt.yticks(sorted_x_values, ytickz, fontsize = Fontsize)
plt.xticks(fontsize=Fontsize)
# plt.axhline(y=np.mean([item[1][0] for item in sorted_acc]), color='gray', linestyle='dotted')
# plt.ylabel("gene sets")
plt.xlabel("AUC change", fontsize = Fontsize)
import matplotlib.ticker as ticker
plt.gca().xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
maxval = np.max(np.abs(sorted_y_values))
plt.xlim(-maxval*1.1,maxval*1.1)
plt.axvline(color = 'grey', linestyle = '--')
plt.title("impact of gene set masking on AUC", fontsize = 12)
# plt.savefig('KEGG_AUC_drop.svg', transparent=True, facecolor='none', format = 'svg')
plt.show()


#### panel b ####
# Load first dataset
with open('data/r_auc_ngenes_Final.pkl','rb') as f:
    r_auc = pkl.load(f)

# Load second dataset
with open('data/r_auc_ngenes_unconnected_Final.pkl','rb') as f:
    r_auc_unconnected = pkl.load(f)

with open('data/n_auc_ngenes_Final.pkl','rb') as f:
    n_auc = pkl.load(f)

# Compute means
r_means = {k: np.mean(v) for k, v in r_auc.items()}
r_means_unconnected = {k: np.mean(v) for k, v in r_auc_unconnected.items()}
n_means = {k: np.mean(v) for k, v in n_auc.items()}

# Compute standard deviations for error bars
r_std = {k: np.std(v) for k, v in r_auc.items()}
r_std_unconnected = {k: np.std(v) for k, v in r_auc_unconnected.items()}

# Sort by number of removed genes
sorted_keys = sorted(r_means.keys())

auc_connected = [r_means[k] for k in sorted_keys]
auc_unconnected = [r_means_unconnected[k] for k in sorted_keys]
auc_n_connected = [n_means[k] for k in sorted_keys]

std_connected = [r_std[k] for k in sorted_keys]
std_unconnected = [r_std_unconnected[k] for k in sorted_keys]


# First data points (baseline values)
baseline_connected = auc_connected[0]
baseline_unconnected = auc_unconnected[0]

# ---- Plot with error bars ----
import matplotlib.ticker as ticker
plt.rcParams['font.family'] = 'Arial'
plt.figure(figsize=(6, 6))


plt.errorbar(
    sorted_keys, auc_connected,
    # yerr=std_connected,
    marker='o',
    capsize=4,
    label='Random gene removal',
    color = "#5E4FA2"
)

# plt.errorbar(
#     sorted_keys, auc_unconnected,
#     yerr=std_unconnected,
#     marker='s',
#     capsize=4,
#     label='Unconnected graph'
# )

plt.errorbar(
    sorted_keys, auc_n_connected,
    marker='s',
    capsize=4,
    label='Sorted gene removal',
    color = "#9E0142"
)
# Dashed horizontal baselines
# plt.axhline(baseline_connected, linestyle='--', linewidth=1, c = "#5E4FA2")
# plt.axhline(baseline_unconnected, linestyle='--', linewidth=1, c = "#9E0142")
ax = plt.gca()
ax.set_ylim(top=1.0)  # ensure 1.0 is included
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))

plt.xlabel("Removed genes")
plt.ylabel("ROC AUC")
plt.title("Impact of gene masking on ROC AUC")
plt.legend()
plt.tight_layout()
# plt.savefig('E_gene_removal.svg', transparent=True, facecolor='none', format = 'svg')
plt.show()