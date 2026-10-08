"""Prepare TCGA gene-expression samples for modeling.

This script loads normalized expression data, filters to suitable cancer projects,
creates a binary label for mutation/CNV-positive samples, and converts each
project into graph-like structures for downstream spatial GCNN SPOT-Mut training/testing/validation.
"""

import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import RobustScaler
import torch

# Read normalized expression matrix: rows = genes, columns = samples.
# Transpose so that rows become samples and columns become genes.
dataset = pd.read_csv("data/LogNormalized_ALL_counts.csv", index_col=0).T

# Read TCGA clinical/annotation metadata.
metadata = pd.read_csv("data/ALL_TCGA_meta.csv", index_col=0)

# Restrict metadata to the same sample IDs appearing in the expression matrix.
metadata = metadata.loc[dataset.index]

# Keep only projects with enough samples and a reasonable mutation rate.
# This avoids highly imbalanced or tiny cohorts that would be poor training targets.
projects = pd.factorize(metadata['cases_project.project_id'])[1]
pass_projects = {}
for Proj in projects:
    tmp_meta = metadata[metadata['cases_project.project_id'] == Proj]
    n_samples = len(tmp_meta.index)
    n_mutated = sum(tmp_meta.loc[:, 'mutated'])
    if n_mutated / n_samples >= 0.1 and n_mutated / n_samples <= 0.85 and n_samples > 250:
        pass_projects[Proj] = (n_samples, n_mutated, n_samples - n_mutated, n_mutated / n_samples)
        print(f"project: {Proj} \t\t #samples: {n_samples}\t #mutated: {n_mutated} ({(n_mutated / n_samples) * 100:.2f}%)\t #WT: {n_samples - n_mutated}")

# Filter the dataset and metadata to the accepted projects only.
metadata = metadata[metadata['cases_project.project_id'].isin(pass_projects.keys())]
dataset = dataset.loc[metadata.index]

# Define the binary target as mutated OR CNV loss.
dataset['Y'] = (metadata['mutated'] | metadata['CNV_loss']).astype(int)

# Split by patient ID so that samples from the same patient do not leak across train/test.
patients = metadata['cases_submitter_id'].unique()
train_samples, test_samples = train_test_split(patients, train_size=0.7, random_state=45)  # 12913

meta_train = metadata[metadata['cases_submitter_id'].isin(train_samples)]
meta_test = metadata[metadata['cases_submitter_id'].isin(test_samples)]
training_data = dataset.loc[meta_train.index]
testing_data = dataset.loc[meta_test.index]

# Scale gene expression features using a robust scaler fitted on the training set only.
# This helps reduce the effect of outliers while preserving relative expression structure.
Y_train = training_data.pop('Y')
scaler = RobustScaler().fit(training_data)
training_data = pd.DataFrame(scaler.transform(training_data), columns=dataset.columns[:-1], index=Y_train.index)
training_data = pd.concat([training_data, Y_train], axis=1)

Y_test = testing_data.pop('Y')
testing_data = pd.DataFrame(scaler.transform(testing_data), columns=dataset.columns[:-1], index=Y_test.index)
testing_data = pd.concat([testing_data, Y_test], axis=1)

# Convert each project into spatial-graph training/test inputs.
# The downstream model expects graph-structured tensor data produced by CreateSTGraph.
from functions.CreateSTGraph import CreateSTGraph

test_data = {}
data = {}
for Project in pass_projects.keys():
    print(Project)
    data[Project] = CreateSTGraph(
        training_data[meta_train['cases_project.project_id'] == Project],
        grid_size=65,
        num_circles=7,
        circle_radius=[6, 8],
        seed=120111,
    )[0]
    test_data[Project] = CreateSTGraph(
        testing_data[meta_test['cases_project.project_id'] == Project],
        grid_size=65,
        num_circles=7,
        circle_radius=[5, 12],
    )[0]
