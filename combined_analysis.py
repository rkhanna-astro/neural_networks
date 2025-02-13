import network_analysis
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker

# Importing fits file
file_name = 'CEREBELLUM40.fits'
cerebellum_properties = network_analysis.network_properties(file_name)

file_name = "map_gasDMD_188_1024_X_D.fits"
cosmic_web_properties = network_analysis.network_properties(file_name)

file_name = 'CORTEX40.fits'
cortex_properties = network_analysis.network_properties(file_name)

cer_total_nodes = cerebellum_properties['total_nodes']
cer_clustering_coefficient = cerebellum_properties['clustering_coefficient']
cer_degree_centrality = cerebellum_properties['degree_centrality']
cer_connectivity = cerebellum_properties['average_connections']

cor_total_nodes = cortex_properties['total_nodes']
cor_clustering_coefficient = cortex_properties['clustering_coefficient']
cor_degree_centrality = cortex_properties['degree_centrality']
cor_connectivity = cortex_properties['average_connections']

cos_total_nodes = cosmic_web_properties['total_nodes']
cos_clustering_coefficient = cosmic_web_properties['clustering_coefficient']
cos_degree_centrality = cosmic_web_properties['degree_centrality']
cos_connectivity = cosmic_web_properties['average_connections']

# This code is just for plotting the histograms as done in figure 3
fig, axes = plt.subplots(2)

cer_counts, cer_bins = np.histogram(cer_clustering_coefficient, bins=10)
cer_fractions_clustering = cer_counts / cer_total_nodes

cor_counts, cor_bins = np.histogram(cor_clustering_coefficient, bins=10)
cor_fractions_clustering = cor_counts / cor_total_nodes

cos_counts, cos_bins = np.histogram(cos_clustering_coefficient, bins=10)
cos_fractions_clustering = cos_counts / cos_total_nodes

axes[0].set_title("Clustering Coefficient")
axes[0].set_xlabel("clustering coefficient C")
axes[0].set_ylabel(r'N(C)/$N_{\text{total}}$')
axes[0].hist(cer_bins[:-1], cer_bins, weights=cer_fractions_clustering, histtype='step', color= 'red', linewidth = 2, label='Cerebellum')
axes[0].hist(cor_bins[:-1], cor_bins, weights=cor_fractions_clustering, histtype='step', color= 'orange', linewidth = 2, label='Cortex')
axes[0].hist(cos_bins[:-1], cos_bins, weights=cos_fractions_clustering, histtype='step', color= 'blue', linewidth = 2, label='Cosmic Web')
axes[0].set_xlim(0, 0.8)
axes[0].set_ylim(0, 0.6)
axes[0].legend(loc='upper right')
axes[0].minorticks_on()

axes[0].xaxis.set_minor_locator(ticker.MultipleLocator(0.1))  # Minor ticks every 0.5 units
axes[0].yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
axes[0].grid(which='both', linestyle='--', linewidth=0.5)

cer_counts_2, cer_bins_2 = np.histogram(cer_degree_centrality, bins=10)
cer_fractions_degree = cer_counts_2 / cer_total_nodes

cor_counts_2, cor_bins_2 = np.histogram(cor_degree_centrality, bins=10)
cor_fractions_degree = cor_counts_2 / cor_total_nodes

cos_counts_2, cos_bins_2 = np.histogram(cos_degree_centrality, bins=10)
cos_fractions_degree = cos_counts_2 / cos_total_nodes

axes[1].set_title("Degree Centrality")
axes[1].set_xlabel(r'degree centrality $C_d$')
axes[1].set_ylabel(r'$N(C_d)$/$N_{\text{total}}$')
axes[1].hist(cer_bins_2[:-1], cer_bins_2, weights=cer_fractions_degree, color= 'red', histtype='step', linewidth = 2, label='Cerebellum')
axes[1].hist(cor_bins_2[:-1], cor_bins_2, weights=cor_fractions_degree, color= 'orange', histtype='step', linewidth = 2, label='Cortex')
axes[1].hist(cos_bins_2[:-1], cos_bins_2, weights=cos_fractions_degree, color= 'blue', histtype='step', linewidth = 2, label='Cosmic Web')
axes[1].set_xlim(0, 0.006) 
axes[1].set_ylim(0, 0.6)
axes[1].legend(loc='upper right')

axes[1].minorticks_on()

axes[1].xaxis.set_minor_locator(ticker.MultipleLocator(0.001))  # Minor ticks every 0.5 units
axes[1].yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
axes[1].grid(which='both', linestyle='--', linewidth=0.5)

# plt.hist(clustering_coefficient, bins=30, density=True, color='skyblue', edgecolor='black', alpha=0.7)
# plt.imshow(dark_matter_density_data, cmap='gray', origin='lower')
# plt.colorbar()
# plt.title("FITS Image")

plt.tight_layout()
plt.show()
