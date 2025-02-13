## From Majd
## This page is for Midterm project of Neural network.
'''
It is based on the paper "The Quantitative comparison between the Neuronal Network and the Cosmic Web" by
"F. Vazza and A. Feletti". A link to the paper
is https://www.frontiersin.org/journals/physics/articles/10.3389/fphy.2020.525731/full. #All the data used in the
paper are in the website for the authors
'''


#This is a python code to mimic a figure in the paper. The figure that we are aimin to mimic is figure 3 from the paper
## Importing the libraries and modules

import matplotlib.image as mpimg
from IPython.display import Image
import numpy as np
#import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import stats
from scipy.stats import shapiro, kstest
from scipy.ndimage import label
from scipy.ndimage import sum_labels


## Themimage we want to mimic is the following

# Image(filename="Figure_to_mimic_j.jpg", width=500, height=500) # the image to mimic is called "Figure_to_mimic_j.jpg"
# img = mpimg.imread("Figure_to_mimic_j.jpg")

# plt.figure(1)
# plt.imshow(img)
# plt.axis("off")  # Hide axes
# plt.show()

## The code for the fits file

def extract_nodes(file_name: str, threshold):

    if 'map' in file_name:
        # the file name. the file should be in the same directory
        fits_file = fits.open(file_name) # opening the fits file
        image_data = fits_file[0].data[1]

    else:
        fits_file = fits.open(file_name) # opening the fits file
        image_data = fits_file[0].data
    
    max_density = np.max(image_data)
    length = image_data.shape[0]
    width = image_data.shape[1]

    # Normalizing the values for lower computational costs
    # for x in range(length):
    #     for y in range(width):
    #         image_data[x][y] = image_data[x][y]/max_density

    # image_data = image_data[1] # 0 for baryonic matter, 1 for Dark matter web
    # This is the threshold for the nodes
    # Creating a binary array of 1 and 0 (True and False). Whenever we have a pixel with a value of more
    # than the threshold, this will create 1 (True), and whenever the pixel is below the threshold, it will
    # create a value of 0

    binary_mask = image_data > threshold
    # We use the "label" function from scipy library to classify the pixels together as one node.
    # The node here represents a cluster of pixels, not necessarily a pixel by itself.

    labeled_array, num_nodes = label(binary_mask)
    '''labeled_array is a matrix to distinguish the pixels
    that are together from the pixels that are not together, based on the location of the pixels.
    If the pixel has a value above the threshold put, and the pixel on its right/left/below/above also have
    a value that is higher than the threshold, then it is considered as one Node. This process keeps
    going till the pixel that is on the top, right, left and (and not or) bottom, have a value less than the threshold 
    num_modes is a values that prints the number of nodes!
    '''
    # Summing the pixel value in each node
    node_fluxes = sum_labels(image_data, labeled_array, index=np.arange(1, num_nodes + 1))
    ## Extracting the X and Y coordinates for each node
    nodes = {}  # Dictionary to store node information

    for node_id in range(1, num_nodes + 1):  # Skip background (0)
        y_coords, x_coords = np.where(labeled_array == node_id)  # Find pixel locations for each node
        total_flux = node_fluxes[node_id - 1]  # Get total flux for this node, by combining all the pixels of the node
        nodes[node_id] = {
            "coords": list(zip(x_coords, y_coords)),  # Store coordinates
            "total_flux": total_flux
        }
    # number_of_nodes_to_print = 3 # how many nodes to print

    # for node_id, data in list(nodes.items())[:number_of_nodes_to_print]:  # Print first 5 nodes
    #     print(f"Node {node_id}:")
    #     print(f"  Coordinates: {data['coords']}")
    #     print(f"  Total Flux: {data['total_flux']:.3e}\n")
    print(f"the total number of nodes with a value higher than threshold {threshold} is " + str(num_nodes))

    # Plot the labeled nodes
    plt.figure(2)
    plt.imshow(labeled_array, cmap='nipy_spectral', origin='lower')
    plt.colorbar(label="Node Label")
    plt.title(f"Identified Nodes (Threshold: {threshold})")
    plt.show()

    return {
        'nodes': nodes, 
        'plots': labeled_array,
    }
