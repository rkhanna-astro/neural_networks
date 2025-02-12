# From pratik
# import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from astropy.io import fits
from astropy import table

from scipy.fft import fft2, fftshift
import scipy.ndimage
# import networkx as nx
from scipy.ndimage import gaussian_filter

## Read in the data
brain_data = fits.getdata('CEREBELLUM40.fits').astype('float64')

smoothed_data = gaussian_filter(brain_data, sigma=2)
downsample_data = smoothed_data[::2, ::2]

print(f'Orginal Data shape: {brain_data.shape}\nDownsampled Data shape: {downsample_data.shape}')

plt.imshow(downsample_data)
plt.show()

# data = downsample_data

data = brain_data

# Step 1: Identify high-intensity peaks (top 10% pixels)
thresh = np.percentile(data, 99)
peaks = np.argwhere(data >= thresh)

# Step 2: Compute enclosed average intensity until threshold delta is met
def find_radius(data, peak, delta):
    x, y = peak
    max_radius = 48 # 16 micrometers = 48 pixel after downsampling
    # max_radius = min(data.shape) // 80  # Limit search radius
    for r in range(1, max_radius):
        mask = (np.arange(data.shape[0])[:, None] - x) ** 2 + (np.arange(data.shape[1])[None, :] - y) ** 2 <= r ** 2
        avg_intensity = np.mean(data[mask])
        if avg_intensity <= delta:
            return r
    return max_radius  # If delta is never met, assign max radius

# Define delta threshold
Delta = np.percentile(data, 98)  # intensity as threshold

# Compute node radii
nodes = []
radii = []
for peak in peaks:
    r = find_radius(data, peak, Delta)
    nodes.append((peak[0], peak[1]))
    radii.append(r)
    
# Step 3: Construct adjacency matrix
num_nodes = len(nodes)
adjacency_matrix = np.zeros((num_nodes, num_nodes), dtype='int32')

for i in range(num_nodes):
    for j in range(i + 1, num_nodes):
        dist = np.linalg.norm(np.array(nodes[i]) - np.array(nodes[j]))
        if dist <= (radii[i] + radii[j]):
            adjacency_matrix[i, j] = 1
            adjacency_matrix[j, i] = 1
            
            
node_degrees = np.sum(adjacency_matrix, axis=1)
clustering_coeffs = []

def compute_clustering_coefficient(adj_matrix, node_idx):
    neighbors = np.where(adj_matrix[node_idx] == 1)[0]
    if len(neighbors) < 2:
        return 0.0
    possible_links = len(neighbors) * (len(neighbors) - 1) / 2
    actual_links = sum(adj_matrix[i, j] for i in neighbors for j in neighbors if i < j)
    return actual_links / possible_links if possible_links > 0 else 0.0

for i in range(num_nodes):
    clustering_coeffs.append(compute_clustering_coefficient(adjacency_matrix, i))

print("Number of nodes:", num_nodes)
print("Number of edges:", np.sum(adjacency_matrix) // 2)
print("Average degree:", np.mean(node_degrees))
print("Average clustering coefficient:", np.mean(clustering_coeffs))

# Step 5: Plot distributions
fig, ax = plt.subplots(1, 2, figsize=(12, 5))
ax[0].hist(clustering_coeffs, bins=20, color='b', alpha=0.7)
ax[0].set_title("Clustering Coefficient Distribution")
ax[0].set_xlabel("Clustering Coefficient")
ax[0].set_ylabel("Frequency")

ax[1].hist(node_degrees, bins=20, color='r', alpha=0.7)
ax[1].set_title("Degree Distribution")
ax[1].set_xlabel("Degree")
ax[1].set_ylabel("Frequency")

plt.tight_layout()
plt.show()

print("The reconstructed connections do not take into account the long-range neural connections, and the clusters shown are purely spatial.")
