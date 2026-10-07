import numpy as np
import torch
from torch_geometric.data import Data
import matplotlib.pyplot as plt
def CreateSTGraph(ST, num_circles = 4, grid_size = 65, circle_radius = [4,8], ones = False, n_slides = 1, seed = np.nan):
    ## Function to create a new spatial transcriptomics (ST) slide with each spot being represented as a node in a graph.
    ## Each node in the Graph is connected to it's 8 neighbouring node (or less if it is at the edge of the ST slide).
    ## Each node contains a transcriptional profile taken from the 'dataset' input variable.
    ## This 'dataset' variable should be a normalized expression matrix loaded as a pandas dataframe with the rows being the samples and the columns being the genes
    ## Last column should be the binary target
    ## precompute datasets
    data1 = ST[ST.iloc[:, -1] == 1]
    data1_idx = len(data1.index)
    data0 = ST[ST.iloc[:, -1] == 0]
    data0_idx = len(data0.index)
    data = {}
    grid = grid_size*grid_size
    if not np.isnan(seed):
        np.random.seed(seed)
        
    for n_slide in range(n_slides):
        samples = []
        # Create node features and labels
        node_features = torch.zeros(grid, len(ST.columns)-1)  # Node features
        node_labels = torch.zeros(grid, dtype=torch.float)  # Node labels
        radii = np.random.choice(range(circle_radius[0], circle_radius[1]), size=num_circles)
        
        # Define circle midpoints
        circle_midpoints = []
        for _ in range(num_circles):
            midpoint = (np.random.choice(range(circle_radius[0], grid_size - circle_radius[0]+1)), np.random.choice(range(circle_radius[0], grid_size - circle_radius[0]+1)))
            circle_midpoints.append(midpoint)

        # Random scramble the data
        index_1s = np.random.choice(data1_idx, replace = False, size = data1_idx)
        index_1s_ind = 0
        index_0s = np.random.choice(data0_idx, replace = False, size = data0_idx)
        index_0s_ind = 0
        new_seeds = np.random.choice(grid, replace = False, size = int(np.ceil(((grid)/data1_idx)+((grid)/data0_idx))))
        new_seeds_ind = 0
        
        # Assign class labels based on circle definitions
        for i in range(grid_size):
            for j in range(grid_size):
                node_idx = i * grid_size + j
                x = torch.tensor([i, j])

                # Check if the node is inside any circle
                inside_circle = False
                for d,midpoint in enumerate(circle_midpoints):
                    center = torch.tensor(midpoint, dtype=torch.float)
                    distance = torch.norm(x.float() - center)
                    if distance <= radii[d]:
                        inside_circle = True
                        break

                if inside_circle or ones:
                    node_labels[node_idx] = 1
                    # Sample features from the distribution for class 1 nodes
                    random_index = index_1s[index_1s_ind]
                    index_1s_ind += 1
                    samples.append(data1.index[random_index])
                    node_features[node_idx] = torch.tensor(data1.iloc[random_index,:-1].values)
                    
                else:
                    # Sample features from the distribution for class 0 nodes
                    random_index = index_0s[index_0s_ind]
                    index_0s_ind += 1
                    samples.append(data0.index[random_index])
                    node_features[node_idx] = torch.tensor(data0.iloc[random_index,:-1].values)
                
                # regenerate new random list on indices
                if index_0s_ind == data0_idx: 
                    if not np.isnan(seed):
                        np.random.seed(new_seeds[new_seeds_ind])
                        new_seeds_ind += 1
                    index_0s = np.random.choice(data0_idx, replace = False, size = data0_idx)
                    index_0s_ind = 0
                if index_1s_ind == data1_idx: 
                    if not np.isnan(seed):
                        np.random.seed(new_seeds[new_seeds_ind])
                        new_seeds_ind += 1
                    index_1s = np.random.choice(data1_idx, replace = False, size = data1_idx)
                    index_1s_ind = 0


        # Create edge index
        edge_index = []
        for i in range(grid_size):
            for j in range(grid_size):
                node_idx = i * grid_size + j
                if i > 0:
                    edge_index.append([node_idx, (i - 1) * grid_size + j])  # Upper neighbor
                if i < grid_size - 1:
                    edge_index.append([node_idx, (i + 1) * grid_size + j])  # Lower neighbor
                if j > 0:
                    edge_index.append([node_idx, i * grid_size + (j - 1)])  # Left neighbor
                if j < grid_size - 1:
                    edge_index.append([node_idx, i * grid_size + (j + 1)])  # Right neighbor
                    # Add edges with diagonal nodes
                if i > 0 and j > 0:
                    edge_index.append([node_idx, (i - 1) * grid_size + (j - 1)])
                if i > 0 and j < grid_size - 1:
                    edge_index.append([node_idx, (i - 1) * grid_size + (j + 1)])
                if i < grid_size - 1 and j > 0:
                    edge_index.append([node_idx, (i + 1) * grid_size + (j - 1)])
                if i < grid_size - 1 and j < grid_size - 1:
                    edge_index.append([node_idx, (i + 1) * grid_size + (j + 1)])

                

        edge_index = torch.tensor(edge_index).t().contiguous()

        # Create the Data object
        data[n_slide] = Data(x=node_features, edge_index=edge_index, y=node_labels, samples = samples)#torch.ones(64*64, dtype=torch.bool))#
    return data