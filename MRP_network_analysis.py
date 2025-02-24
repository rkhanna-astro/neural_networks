'''
It is based on the paper "The Quantitative comparison between the Neuronal Network and the Cosmic Web" by
"F. Vazza and A. Feletti". A link to the paper
is https://www.frontiersin.org/journals/physics/articles/10.3389/fphy.2020.525731/full. #All the data used in the
paper are in the website for the authors
'''

import numpy as np
from astropy.io import fits
import select_and_cluster_pixels as extract_nodes

'''
This function calculates the degree centrality for each node (or pixel).
For a particular node, it is calculated by calculating the number of nodes it is connected to
and dividng it by total nodes - 1 (itself).  
'''
def calculate_degree_centrality(adjacency_matrix):
    total_nodes = len(adjacency_matrix[0])
    degree_centrality = [0]*total_nodes

    for x in range(total_nodes):
        edges = sum(adjacency_matrix[x])
        degree_centrality[x] = edges/(total_nodes-1)

    return degree_centrality

'''
This function calculates the clustering coefficient for each node (or pixel).
For a particular node, we first calculate the number of triangles it forms with the pair of nodes 
its connected to. Then we divide that number with the maximum possible triangle that the node can be part
of with its connected nodes. (nC2 or n(n-1)/2)  
'''
def calculate_clustering_coefficient(adjacency_matrix):
    total_nodes = len(adjacency_matrix[0])
    clustering_coefficient = [0]*total_nodes

    '''
    Here we compute the clustering coefficient for each node
    This calculates ratio of the no. of triangles that our each node
    form compared to the total triangles it can form based on the edges it has.
    '''
    for x in range(total_nodes):
        edges = sum(adjacency_matrix[x])

        if edges < 2:
            continue

        traingles = 0
        for y in range(total_nodes):
            if adjacency_matrix[x][y] == 0:
                continue
            
            for z in range(y+1, total_nodes):
                if adjacency_matrix[x][z] == 1 and adjacency_matrix[y][z] == 1:
                    traingles += 1

        clustering_coefficient[x] = 2*traingles/(edges*(edges-1))
    
    return clustering_coefficient

'''
This function calculates the total edges that our network has.
'''
def calculate_total_edges(adjacency_matrix):
    total_nodes = len(adjacency_matrix[0])
    total_edges = 0

    for x in range(total_nodes):
        edges = sum(adjacency_matrix[x])
        total_edges += edges
    
    return total_edges

'''
A particular node from the dictionary returned by select_and_cluster_pixel module has a collection of pixels
clustered together that are above threshold in intensity value. But we only one pixel to represent that
collection of clustered pixels, for further creation of edges.

The function below has the logic to find the centroid pixel from the set of clustered pixels to represent
our set of selected nodes. In the end, all these centroid nodes will be used to form edges.
'''
def find_potential_nodes(nodes):
    potential_nodes = set()

    for node, values in nodes.items():
        coordinates = values['coords']
        sum_x = 0
        sum_y = 0

        for x, y in coordinates:
            sum_x += x
            sum_y += y

        avg_x = sum_x/len(coordinates)
        avg_y = sum_y/len(coordinates)

        initial_distance = float('inf')

        cen_x = 0
        cen_y = 0

        for x, y in coordinates:
            distance = (avg_x - x)**2 + (avg_y - y)**2

            if distance < initial_distance:
                initial_distance = distance
                cen_x = x
                cen_y = y

        # We add these nodes to the potential nodes set for then to create adjacency matrix.
        potential_nodes.add((cen_x, cen_y))
    
    return potential_nodes

'''
This is an another way of selecting pixels above threshold. Here we first take all the pixels above threshold.
We then select a random pixel to be part of final nodes list and remove all the neighbouring nodes 
using breadth-first search algorithm. We didn't use this algorithm but it was giving comparable results 
when compared to SciPy clustering algorithm.
'''
def breadth_first_selection_nodes(pixel_data, threshold):
    length = len(pixel_data)
    width = len(pixel_data[0])
    potential_nodes = set()

    for x in range(length):
            for y in range(width):
                if pixel_data[x][y] >= threshold:
                    potential_nodes.add((x,y))

    '''
    After we get a set potential nodes with intensities higher than threshold
    the algorithm removes all the neighbours of a potential node till it can find any.
    These neighbors should also be part of the potential nodes set otherwise it will skip removing.
    It continues till it can find all the neighbors of a initial potential node.
    '''
    for x, y in potential_nodes.copy():
        stack = [(x-1, y), (x+1, y), (x, y+1), (x, y-1)]
        while stack:
            i, j = stack.pop()
            if (i, j) != (x, y) and (i, j) in potential_nodes:
                potential_nodes.discard((i, j))
                stack.append((i + 1, j))
                stack.append((i, j + 1))
                stack.append((i - 1, j))
                stack.append((i, j - 1))
    
    return potential_nodes

'''
The function calculates the network properties for the three types of network studied i.e.
cosmic web, cerebellum and cortex.
'''
def network_properties(file_name: str):
    threshold = None
    flattened_data = None
    cut_off = None
    linking_length = None

    # Importing fits file
    if 'gas' in file_name:
        cosmic_web_intensity_data = fits.open(file_name)
    
        # There are two levels in cosmic data, 1st gas density and 2nd dark matter density.
        # We don't use gas density because of very low intensity compared to dark matter which is much more in amount.
        dark_matter_density_data = cosmic_web_intensity_data[0].data[1]

        # Flattening cosmic web data to get the top 1 % percentile values
        flattened_data = dark_matter_density_data.flatten()
    elif 'CORTEX' in file_name:
        cortex_data = fits.open(file_name)[0].data

        # Flattening cortex data to get the top 0.2 % percentile values
        flattened_data = cortex_data.flatten()
    else:
        cerebellum_data = fits.open(file_name)[0].data
        
        # Flattening cerebellum data to get the top 0.2 % percentile values
        flattened_data = cerebellum_data.flatten()

    '''
    Here, after some hit-and-trial approach we select the threshold percentiles for
    selecting the nodes using SciPy clustering algorithm. Authors in their research
    used a different selecting algorithm called halo-finding. Hence, for our algorithm we had
    to use these threshold values to get same amount of nodes and average connectivity as the
    authors.
    '''
    if 'gas' in file_name:
        cut_off = 98.5
    elif 'CORTEX' in file_name:
        cut_off = 99.8
    else:
        cut_off = 99.8
    
    # print("Here is the cut_off", cut_off)
    threshold = np.percentile(flattened_data, cut_off)
    values = extract_nodes.extract_and_cluster_pixels(file_name, threshold)
    nodes = values['nodes']
    # print("The total number of nodes", len(nodes))

    # Now we can create a set of potential nodes using their co-ordinates in pixel 2-D array.
    selected_nodes = find_potential_nodes(nodes)
    # print(len(selected_nodes))

    # We calculate total number of nodes that we finally select for then to create edges and adjacency matrix.
    number_of_nodes = len(selected_nodes)
    
    # We initiate adjacency square matrix that will be later used to store the edges data.
    adjacency_matrix = np.zeros(shape=(number_of_nodes, number_of_nodes))

    # Now we mark all the selected nodes to a particular position in our matrix
    matrix_map = {}
    pos = 0

    # This is random mapping but will not impact our calculations
    for x, y in selected_nodes:
        matrix_map[(x, y)] = pos
        pos += 1

    '''
    Now we mark the nodes that are close in distance with edges. This is done by controlling a network parameter called the linking length.
    For example, for cosmic web each node is 0.04 Mpc long and we form an edge if the distance between two nodes is less than
    1.2 Mpc. This means we can go in any direction at most 30 co-ordinates (excluding original) in the adjacency matrix.
    This algorithm will check the difference in position of all nodes with current node and if less than 30 co-ordinates
    distance it will form an undirected edge with that node.

    However, linking length was a arbitary parameter as suggested by the authors. For cosmic web and cerebellum we use the
    same linking length values as suggested by authors. But for cortex data, we had to use different value to match the data
    as found in the paper.
    '''
    if 'gas' in file_name:
        linking_length = 30
    elif 'CORTEX' in file_name:
        # This is a bit arbitary
        linking_length = 30
    else:
        linking_length = 15

    for x, y in selected_nodes:
        pos_i = matrix_map[(x, y)]
        for i, j in selected_nodes:
            if i != x and j != y and (abs(i-x)**2 + abs(j-y)**2)**(0.5) < linking_length:
                pos_f = matrix_map[(i, j)]
                adjacency_matrix[pos_i][pos_f] = 1
                adjacency_matrix[pos_f][pos_i] = 1

    '''
    Now we compute the network properites such as degree centrality and clustering coefficient for each node.

    We also compute the general network properties such as total edges and average connection for each node to compare
    our algorithm's output with the properties mentioned in the paper.

    For all the three FITS data files, we were able to achieve very similar general network properties as mentioned in the paper.
    '''
    degree_centrality = calculate_degree_centrality(adjacency_matrix)
    total_edges = calculate_total_edges(adjacency_matrix)
    clustering_coefficient = calculate_clustering_coefficient(adjacency_matrix)
    average_connections = total_edges/number_of_nodes
    # the average value in the paper for cosmic web is 3.8 ~ 4.1

    '''
    From here, we return the network properties to the class or module that imported this module for further analysis or plotting.

    We return the properites in form of dictionary for easier access and better usability or for further extensibility.
    '''
    network_properties = {
        'total_nodes' : number_of_nodes,
        'average_connections': average_connections,
        'degree_centrality': degree_centrality,
        'clustering_coefficient': clustering_coefficient
    }

    return network_properties