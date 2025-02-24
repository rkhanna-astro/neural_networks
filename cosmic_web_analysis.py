'''
It is based on the paper "The Quantitative comparison between the Neuronal Network and the Cosmic Web" by
"F. Vazza and A. Feletti". A link to the paper
is https://www.frontiersin.org/journals/physics/articles/10.3389/fphy.2020.525731/full. #All the data used in the
paper are in the website for the authors
'''

import MRP_network_analysis
import matplotlib.pyplot as plt
import numpy as np

# Importing Cosmic web FITS files for network analysis.
files = ['D']
nodes_data = []
connectivity_data = []

'''
Here, we go through each fits file for COSMIC web images to study their network properties.
We also collect and compare their total nodes data and average connectivity values to see 
how consistent our algorithm is performing for every case.
'''
for file in files:
    file_name = f"map_gasDMD_188_1024_X_{file}.fits"
    cosmic_web_properties = MRP_network_analysis.network_properties(file_name)

    total_nodes = cosmic_web_properties['total_nodes']
    nodes_data.append(total_nodes)

    clustering_coefficient = cosmic_web_properties['clustering_coefficient']
    degree_centrality = cosmic_web_properties['degree_centrality']

    connectivity = cosmic_web_properties['average_connections']
    connectivity_data.append(connectivity)

    # This code is just for plotting the histograms as done in figure 3
    fig, axes = plt.subplots(2)

    fig.suptitle(f"Cosmic Web Network Analysis (File {file}) ", fontsize=14)

    counts, bins = np.histogram(clustering_coefficient, bins=10)
    fractions_clustering = counts / total_nodes

    axes[0].set_title("Clustering Coefficient")
    axes[0].set_xlabel("clustering coefficient C")
    axes[0].set_ylabel("N(C)/N_total")
    axes[0].hist(bins[:-1], bins, weights=fractions_clustering, histtype='step', linewidth = 2)
    # axes[0].bar(bins[:-1], fractions_clustering, width=np.diff(bins), edgecolor='blue', facecolor='none', linewidth=2, color='none')

    counts_2, bins_2 = np.histogram(degree_centrality, bins=10)
    fractions_degree = counts_2 / total_nodes

    axes[1].set_title("Degree Centrality")
    axes[1].set_xlabel("degree centrality C_d")
    axes[1].set_ylabel("N(C_d)/N_total")
    axes[1].hist(bins_2[:-1], bins_2, weights=fractions_degree, histtype='step', linewidth = 2)

plt.tight_layout()
plt.show()
