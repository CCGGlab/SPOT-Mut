import torch
import numpy as np
import pickle as pkl

import matplotlib.pyplot as plt

import torch.nn as nn

from scipy.stats import mannwhitneyu, false_discovery_control

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(device)

bootstrap_results = pkl.load(open("../data/Mutation_Classes_100bootstrap_recall_results.pkl", "rb"))

## Plot Code
#### Fig 2c ####

PLOT_ORDER = ['Truncating','Only Missense','True Silent']
group_names = [name for name in PLOT_ORDER if name in bootstrap_results]
median_data = [bootstrap_results[name] for name in group_names]

fig, ax = plt.subplots(figsize=(1.1+1.2*len(group_names), 4.5))
ax.boxplot(median_data, tick_labels=group_names, showfliers=False, widths = 0.65, medianprops=dict({'color': '#9E0142'}))
for i, d in enumerate(median_data):
    x = np.random.normal(i + 1, 0.08, size=len(d))
    ax.scatter(x, d, alpha=0.75, s=25, color='#5E4FA2', linewidth=0)


ax.set_ylim(0, 1.2)
ax.set_title('Recall for Truncating and Missense Clones')
ax.set_ylabel('Recall')
ax.set_xticklabels(['Truncating','Missense', 'Silent'], rotation=30, ha='right')

# Compute pairwise significance
significance_pairs = {}
for i in range(len(group_names)):
    for j in range(i+1, len(group_names)):
        name1 = group_names[i]
        name2 = group_names[j]
        data1 = bootstrap_results[name1]
        data2 = bootstrap_results[name2]
        _, p = mannwhitneyu(data1, data2, alternative='two-sided')
        significance_pairs[(i, j)] = p

# Apply Bonferroni correction
corrected_ps = false_discovery_control(list(significance_pairs.values()))
corrected_pairs = {item: corrected_ps[i] for i, item in enumerate(significance_pairs.keys())}
print(corrected_pairs)

# Add significance indicators with lines at different heights for corrected p < 0.05
# Define manual order
manual_order = [(0, 1), (1, 2), (0, 2)]
line_height_offset = 0.05
line_spacing = 0.035
cap_height = 0.02
asterisk_offset = -0.01
idx = 0
for pair in manual_order:
    if pair in corrected_pairs and corrected_pairs[pair] < 0.05:
        i, j = pair
        y_line = max(np.max(median_data[i]), np.max(median_data[j])) + line_height_offset + line_spacing * idx
        # Draw horizontal line
        ax.plot([i+1, j+1], [y_line, y_line], color='black', linewidth=1)
        # Draw end caps (small vertical lines)
        ax.plot([i+1, i+1], [y_line-0.001, y_line - cap_height], color='black', linewidth=1)
        ax.plot([j+1, j+1], [y_line-0.001, y_line - cap_height], color='black', linewidth=1)
        # Draw asterisk
        ax.text((i + j + 2) / 2, y_line + asterisk_offset, '*', ha='center', va='bottom', fontsize=14, color='black')
        idx += 1
plt.rcParams['font.family'] = 'Arial'
plt.tight_layout()
# plt.savefig('C_boxplot_trunc_missense.svg', transparent=True, facecolor='none', format = 'svg')
plt.show()

#### Fig 2d ####
PLOT_ORDER = ['Double Hit', 'Single Hit', 'WT CNV LOH']
group_names = [name for name in PLOT_ORDER if name in bootstrap_results]
median_data = [bootstrap_results[name] for name in group_names]

fig, ax = plt.subplots(figsize=(1.1+1.2*len(group_names), 4.5))
ax.boxplot(median_data, tick_labels=group_names, showfliers=False, widths = 0.65, medianprops=dict({'color': '#9E0142'}))
for i, d in enumerate(median_data):
    x = np.random.normal(i + 1, 0.08, size=len(d))
    ax.scatter(x, d, alpha=0.75, s=25, color='#5E4FA2', linewidth=0)


ax.set_ylim(0, 1.2)
ax.set_title('Recall for Single and Double Hit Clones')
ax.set_ylabel('Recall')
ax.set_xticklabels(['Double Hit', 'Single Hit SNV', 'Single Hit CNV'], rotation=30, ha='right')

# Compute pairwise significance
significance_pairs = {}
for i in range(len(group_names)):
    for j in range(i+1, len(group_names)):
        name1 = group_names[i]
        name2 = group_names[j]
        data1 = bootstrap_results[name1]
        data2 = bootstrap_results[name2]
        _, p = mannwhitneyu(data1, data2, alternative='two-sided')
        significance_pairs[(i, j)] = p

# Apply Bonferroni correction
corrected_ps = false_discovery_control(list(significance_pairs.values()))
corrected_pairs = {item: corrected_ps[i] for i, item in enumerate(significance_pairs.keys())}
print(corrected_pairs)

# Add significance indicators with lines at different heights for corrected p < 0.05
# Define manual order
manual_order = [(0, 2),(0, 1), (1, 2) ]
line_height_offset = 0.05
line_spacing = 0.06
cap_height = 0.02
asterisk_offset = -0.01
idx = 0
for pair in manual_order:
    if pair in corrected_pairs and corrected_pairs[pair] < 0.05:
        i, j = pair
        y_line = max(np.max(median_data[i]), np.max(median_data[j])) + line_height_offset + line_spacing * idx
        # Draw horizontal line
        ax.plot([i+1, j+1], [y_line, y_line], color='black', linewidth=1)
        # Draw end caps (small vertical lines)
        ax.plot([i+1, i+1], [y_line-0.001, y_line - cap_height], color='black', linewidth=1)
        ax.plot([j+1, j+1], [y_line-0.001, y_line - cap_height], color='black', linewidth=1)
        # Draw asterisk
        ax.text((i + j + 2) / 2, y_line + asterisk_offset, '*', ha='center', va='bottom', fontsize=14, color='black')
        idx += 1
plt.rcParams['font.family'] = 'Arial'
plt.tight_layout()
# plt.savefig('D_boxplot_single_double.svg', transparent=True, facecolor='none', format = 'svg')
plt.show()


#### Fig 2e ####
PLOT_ORDER = ['High Impact LOH', 'Low Impact LOH', 'High Impact CNN', 'Low Impact CNN']
group_names = [name for name in PLOT_ORDER if name in bootstrap_results]
median_data = [bootstrap_results[name] for name in group_names]

fig, ax = plt.subplots(figsize=(1.1+1.2*len(group_names), 5.5))
ax.boxplot(median_data, tick_labels=group_names, showfliers=False, widths = 0.65, medianprops=dict({'color': '#9E0142'}))
for i, d in enumerate(median_data):
    x = np.random.normal(i + 1, 0.08, size=len(d))
    ax.scatter(x, d, alpha=0.75, s=25, color='#5E4FA2', linewidth=0)


ax.set_ylim(0, 1.4)
ax.set_title('Recall for Low and High Impact Hemizygous/CNN Clones')
ax.set_ylabel('Recall')
ax.set_xticklabels(['High Impact Hemizygous','Low Impact Hemizygous', 'High Impact Neutral', 'Low Impact Neutral'], rotation=30, ha='right')

# Compute pairwise significance
significance_pairs = {}
for i in range(len(group_names)):
    for j in range(i+1, len(group_names)):
        name1 = group_names[i]
        name2 = group_names[j]
        data1 = bootstrap_results[name1]
        data2 = bootstrap_results[name2]
        _, p = mannwhitneyu(data1, data2, alternative='two-sided')
        significance_pairs[(i, j)] = p

# Apply Bonferroni correction
corrected_ps = false_discovery_control(list(significance_pairs.values()))
corrected_pairs = {item: corrected_ps[i] for i, item in enumerate(significance_pairs.keys())}
print(corrected_pairs)

# Add significance indicators with lines at different heights for corrected p < 0.05
# Define manual order
manual_order = [(0, 1),(2, 3), (1, 2), (0, 2), (1, 3), (0, 3) ]
line_height_offset = 0.05
line_spacing = 0.06
cap_height = 0.02
asterisk_offset = -0.01
idx = 0
for pair in manual_order:
    if pair in corrected_pairs and corrected_pairs[pair] < 0.05:
        i, j = pair
        y_line = max(np.max(median_data[i]), np.max(median_data[j])) + line_height_offset + line_spacing * idx
        # Draw horizontal line
        ax.plot([i+1, j+1], [y_line, y_line], color='black', linewidth=1)
        # Draw end caps (small vertical lines)
        ax.plot([i+1, i+1], [y_line-0.001, y_line - cap_height], color='black', linewidth=1)
        ax.plot([j+1, j+1], [y_line-0.001, y_line - cap_height], color='black', linewidth=1)
        # Draw asterisk
        ax.text((i + j + 2) / 2, y_line + asterisk_offset, '*', ha='center', va='bottom', fontsize=14, color='black')
        idx += 1
plt.rcParams['font.family'] = 'Arial'
plt.tight_layout()
# plt.savefig('E_boxplot_low_high_loh_cnn.svg', transparent=True, facecolor='none', format = 'svg')
plt.show()


#### Suppl Fig 1d ####
PLOT_ORDER = ['High Purity','Low Purity']
group_names = [name for name in PLOT_ORDER if name in bootstrap_results]
median_data = [bootstrap_results[name] for name in group_names]

fig, ax = plt.subplots(figsize=(1.1+1.2*len(group_names), 4.5))
ax.boxplot(median_data, tick_labels=group_names, showfliers=False, widths = 0.65, medianprops=dict({'color': '#9E0142'}))
for i, d in enumerate(median_data):
    x = np.random.normal(i + 1, 0.08, size=len(d))
    ax.scatter(x, d, alpha=0.75, s=25, color='#5E4FA2', linewidth=0)


ax.set_ylim(0, 1.1)
ax.set_title('Recall for High and Low Tumor Purity Clones')
ax.set_ylabel('Recall')
ax.set_xticklabels(['High Purity','Low Purity'], rotation=30, ha='right')

# Compute pairwise significance
significance_pairs = {}
for i in range(len(group_names)):
    for j in range(i+1, len(group_names)):
        name1 = group_names[i]
        name2 = group_names[j]
        data1 = bootstrap_results[name1]
        data2 = bootstrap_results[name2]
        _, p = mannwhitneyu(data1, data2, alternative='two-sided')
        significance_pairs[(i, j)] = p

# Apply Bonferroni correction
corrected_ps = false_discovery_control(list(significance_pairs.values()))
corrected_pairs = {item: corrected_ps[i] for i, item in enumerate(significance_pairs.keys())}
print(corrected_pairs)

# Add significance indicators with lines at different heights for corrected p < 0.05
# Define manual order
manual_order = [(0, 1) ]
line_height_offset = 0.05
line_spacing = 0.06
cap_height = 0.02
asterisk_offset = -0.01
idx = 0
for pair in manual_order:
    if pair in corrected_pairs and corrected_pairs[pair] < 0.05:
        i, j = pair
        y_line = max(np.max(median_data[i]), np.max(median_data[j])) + line_height_offset + line_spacing * idx
        # Draw horizontal line
        ax.plot([i+1, j+1], [y_line, y_line], color='black', linewidth=1)
        # Draw end caps (small vertical lines)
        ax.plot([i+1, i+1], [y_line-0.001, y_line - cap_height], color='black', linewidth=1)
        ax.plot([j+1, j+1], [y_line-0.001, y_line - cap_height], color='black', linewidth=1)
        # Draw asterisk
        ax.text((i + j + 2) / 2, y_line + asterisk_offset, '*', ha='center', va='bottom', fontsize=14, color='black')
        idx += 1
plt.rcParams['font.family'] = 'Arial'
plt.tight_layout()
# plt.savefig('B_boxplot_high_low_purity.svg', transparent=True, facecolor='none', format = 'svg')
plt.show()

