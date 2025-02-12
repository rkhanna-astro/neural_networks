import network_analysis
import matplotlib.pyplot as plt
import numpy as np

# Importing fits file
files = ['A', 'B', 'C', 'D']
nodes_data = []
connectivity_data = []

fig_1, axes_1 = plt.subplots(2)

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

    counts, bins = np.histogram(clustering_coefficient, bins=10)
    fractions_clustering = counts / total_nodes
    axes[0].bar(bins[:-1], fractions_clustering, width=np.diff(bins), color='skyblue', edgecolor='black', alpha=0.7)

    counts_2, bins_2 = np.histogram(degree_centrality, bins=10)
    fractions_degree = counts_2 / total_nodes
    axes[1].bar(bins_2[:-1], fractions_degree, width=np.diff(bins_2), color='blue', edgecolor='black', alpha=0.7)

    # plt.hist(clustering_coefficient, bins=30, density=True, color='skyblue', edgecolor='black', alpha=0.7)
    # plt.imshow(dark_matter_density_data, cmap='gray', origin='lower')
    # plt.colorbar()
    # plt.title("FITS Image")

axes_1[0].plot(files, nodes_data, marker='s')
axes_1[0].set_title(f"Total Nodes")

axes_1[1].plot(files, connectivity_data, marker='o')
axes_1[1].set_title(f"Average Connectivity")
plt.show()
