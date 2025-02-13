import numpy as np
from astropy.io import fits
import AM9642_Neural_netowk_midterm_project______ as extract_nodes

# Importing fits file
# The function network_properties calculate all the properties for both cosmic 
def network_properties(file_name: str):
    threshold = None
    flattened_data = None
    cut_off = None

    if 'gas' in file_name:
        cosmic_web_intensity_data = fits.open(file_name)
    
        # There are two levels in cosmic data, 1st gas density and 2nd dark matter density
        # gas_density_data = cosmic_web_intensity_data[0].data[0]
        dark_matter_density_data = cosmic_web_intensity_data[0].data[1]

        # Flattening cosmic web data to get the top 1 % percentile values
        flattened_data = dark_matter_density_data.flatten()
    elif 'CORTEX' in file_name:
        # Importing the cortex fits file data
        cortex_data = fits.open(file_name)[0].data

        # Flattening cortex data to get the top 0.2 % percentile values
        flattened_data = cortex_data.flatten()
    else:
        # Importing the cerebellum fits file data
        cerebellum_data = fits.open(file_name)[0].data
        
        # Flattening cerebellum data to get the top 0.2 % percentile values
        flattened_data = cerebellum_data.flatten()

    # This gets the top 1% threshold value
    if 'gas' in file_name:
        cut_off = 98.5
    elif 'CORTEX' in file_name:
        cut_off = 99.8
    else:
        cut_off = 99.8
    
    print("Here is the cut_off", cut_off)
    threshold = np.percentile(flattened_data, cut_off)
    values = extract_nodes.extract_nodes(file_name, threshold)
    nodes = values['nodes']
    # labeled_array = values['plots']
    print("The total number of nodes", len(nodes))

    # Now we can create a set of potential nodes using their co-ordinates in pixel 2-D array.
    potential_nodes = set()

    # This algorithm gets the centroid pixel for each node collection.
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

    # for x in range(length):
    #     for y in range(width):
    #         if dark_matter_density_data[x][y] >= threshold:
    #             potential_nodes.add((x,y))

    # # After we get a set potential nodes with intensities higher than threshold
    # # The algorithm removes all the neighbours of a potential node till it can find any
    # # These neighbors should also be part of the potential nodes set otherwise it will skip removing.
    # # It continues till it can find all the neighbors of a initial potential node.
    # for x, y in potential_nodes.copy():
    #     stack = [(x-1, y), (x+1, y), (x, y+1), (x, y-1)]
    #     while stack:
    #         i, j = stack.pop()
    #         if (i, j) != (x, y) and (i, j) in potential_nodes:
    #             potential_nodes.discard((i, j))
    #             stack.append((i + 1, j))
    #             stack.append((i, j + 1))
    #             stack.append((i - 1, j))
    #             stack.append((i, j - 1))

    # # print(len(potential_nodes))

    # We calculate total number of nodes that we finally select for then to create edges and adjacency matrix.
    number_of_nodes = len(potential_nodes)

    # We initiate adjacency matrix
    adjacency_matrix = np.zeros(shape=(number_of_nodes, number_of_nodes))

    # Now we mark all the selected nodes to a particular position in our matrix
    matrix_map = {}
    pos = 0

    # This is random mapping but will not impact our calculations
    for x,y in potential_nodes:
        matrix_map[(x, y)] = pos
        pos += 1

    # Now we mark the close nodes with edges, this is done by using the linking length in paper
    # as each node is 0.04 Mpc long and we form an edge if the distance between two nodes is less than
    # 1.2 Mpc, it means we can go in any direction at most 30 co-ordinates (excluding original)
    # this algorithm will check the difference in position of all nodes with current node
    # and if less than 30 co-ordinates it will form an undirected edge with that node.
    linking_length = None
    if 'gas' in file_name:
        linking_length = 30
    elif 'CORTEX' in file_name:
        # This is a bit arbitary
        linking_length = 30
    else:
        linking_length = 15

    # Based on the linking length now if the pixel is less than that linking length far from
    # our node in consideration we create an edge for those two pixels.
    for x,y in potential_nodes:
        pos_i = matrix_map[(x, y)]
        for i,j in potential_nodes:
            if i != x and j != y and (abs(i-x)**2 + abs(j-y)**2)**(0.5) < linking_length:
                pos_f = matrix_map[(i, j)]
                adjacency_matrix[pos_i][pos_f] = 1
                adjacency_matrix[pos_f][pos_i] = 1

    # Now we compute the degree centrality for each node
    # It's just the ratio of formed edges to total possible edges (total nodes)
    total_edges = 0
    degree_centrality = [0]*number_of_nodes
    for x in range(number_of_nodes):
        edges = sum(adjacency_matrix[x])
        degree_centrality[x] = edges/(number_of_nodes-1)
        total_edges += sum(adjacency_matrix[x])

    # Here we compute the clustering coefficient for each node
    # Inspired from Pratik's code, this calculates ratio of the no. of triangles that our each node
    # form compared to the total triangles it can form based on the edges it has.
    clustering_coefficient = [0]*number_of_nodes
    for x in range(number_of_nodes):
        edges = sum(adjacency_matrix[x])

        if edges < 2:
            continue

        traingles = 0
        for y in range(number_of_nodes):
            if adjacency_matrix[x][y] == 0:
                continue
            
            for z in range(y+1, number_of_nodes):
                if adjacency_matrix[x][z] == 1 and adjacency_matrix[y][z] == 1:
                    traingles += 1

        clustering_coefficient[x] = 2*traingles/(edges*(edges-1))

    # This calculates the average no. of edges or connectivity for each node
    average_connections = total_edges/number_of_nodes
    # the average value in the paper for cosmic web is 3.8 ~ 4.1
    print(average_connections)

    # From this module, we return all the important network properties.
    network_properties = {
        'total_nodes' : number_of_nodes,
        'average_connections': average_connections,
        'degree_centrality': degree_centrality,
        'clustering_coefficient': clustering_coefficient
    }

    return network_properties