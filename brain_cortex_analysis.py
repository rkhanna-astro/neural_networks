'''
It is based on the paper "The Quantitative comparison between the Neuronal Network and the Cosmic Web" by
"F. Vazza and A. Feletti". A link to the paper
is https://www.frontiersin.org/journals/physics/articles/10.3389/fphy.2020.525731/full. #All the data used in the
paper are in the website for the authors
'''

import MRP_network_analysis
import matplotlib.pyplot as plt
import numpy as np

# Importing CORTEX fits file
file_name = 'CORTEX40.fits'

'''
Here, we go through the fits file for the brain CORTEX images to study the network properties.
'''
cortex_properties = MRP_network_analysis.network_properties(file_name)

total_nodes = cortex_properties['total_nodes']
clustering_coefficient = cortex_properties['clustering_coefficient']
degree_centrality = cortex_properties['degree_centrality']
connectivity = cortex_properties['average_connections']

# This code is just for plotting the histograms as done in figure 3
fig, axes = plt.subplots(2)

counts, bins = np.histogram(clustering_coefficient, bins=10)
fractions_clustering = counts / total_nodes

axes[0].set_title("Clustering Coefficient")
axes[0].set_xlabel("clustering coefficient C")
axes[0].set_ylabel("N(C)/N_total")
axes[0].hist(bins[:-1], bins, weights=fractions_clustering, histtype='step', linewidth = 2)

counts_2, bins_2 = np.histogram(degree_centrality, bins=10)
fractions_degree = counts_2 / total_nodes

axes[1].set_title("Degree Centrality")
axes[1].set_xlabel("degree centrality C_d")
axes[1].set_ylabel("N(C_d)/N_total")
axes[1].hist(bins_2[:-1], bins_2, weights=fractions_degree, histtype='step', linewidth = 2)

plt.tight_layout()
plt.show()
