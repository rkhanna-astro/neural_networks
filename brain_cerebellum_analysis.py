import network_analysis
import matplotlib.pyplot as plt
import numpy as np

# Importing fits file
file_name = 'CEREBELLUM40.fits'
cerebellum_properties = network_analysis.network_properties(file_name)

total_nodes = cerebellum_properties['total_nodes']
clustering_coefficient = cerebellum_properties['clustering_coefficient']
degree_centrality = cerebellum_properties['degree_centrality']
connectivity = cerebellum_properties['average_connections']

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
plt.show()
