import network_analysis
import matplotlib.pyplot as plt
import numpy as np

# Importing fits file
files = ['D']
nodes_data = []
connectivity_data = []

# fig_1, axes_1 = plt.subplots(2)

for file in files:
    file_name = f"map_gasDMD_188_1024_X_{file}.fits"
    cosmic_web_properties = network_analysis.network_properties(file_name)

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

    # plt.hist(clustering_coefficient, bins=30, density=True, color='skyblue', edgecolor='black', alpha=0.7)
    # plt.imshow(dark_matter_density_data, cmap='gray', origin='lower')
    # plt.colorbar()
    # plt.title("FITS Image")

# axes_1[0].plot(files, nodes_data, marker='s')
# axes_1[0].set_title(f"Total Nodes")

# axes_1[1].plot(files, connectivity_data, marker='o')
# axes_1[1].set_title(f"Average Connectivity")
plt.tight_layout()
plt.show()
