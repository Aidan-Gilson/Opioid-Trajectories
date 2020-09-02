# DTW Step 1

'''Generates a pairwise distance matrix for all the trajectories for a given distance metrix.
The matrix is later used to cluster the trajectories according to some threshold.'''

import numpy as np
import time
import pickle
import json
from scipy import spatial

# The trajectores found from the greedy algorithm (of all lengths) after convertion to CUI according to the 
# USMLE metathesaurus. 
with open('CUI Trajectories 2.txt') as outfile:
    trajectories = json.load(outfile)

# A file that maps a given CUI code to a 512 dimmensional feature space embedding.
file = open('davidcui2vec.pkl', 'rb')
cui2vec = pickle.load(file)

# Performs the dtw on a pair of trajectories
def dtw(s, t):
    n, m = len(s), len(t)
    # Initializes matrix for the comparison
    dtw_matrix = np.ones((n + 1, m + 1)) * np.inf
    dtw_matrix[0, :] = 0
    dtw_matrix[:, 0] = 0

    # Iterates through the matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            first_vec = cui2vec[s[i - 1]]
            second_vec = cui2vec[t[j - 1]]
            # Changing the cost function changes the distance matrix.
            cost = np.linalg.norm(first_vec-second_vec)
            # cost = spatial.distance.chebyshev(first_vec, second_vec)
            last_min = np.min([dtw_matrix[i - 1, j], dtw_matrix[i, j - 1], dtw_matrix[i - 1, j - 1]])
            dtw_matrix[i, j] = cost + last_min
    # Returns the distance between the two trajectories
    return dtw_matrix[n, m]

# The final matrix that will contain the distance between any two trajectoreis
distance_matrix = np.zeros((len(trajectories), len(trajectories)))

for x in range(len(trajectories)):
    for y in range(x + 1, len(trajectories)):

        # Gets the disntace for those two trajectories and saves them
        distance = dtw(trajectories[x], trajectories[y])
        distance_matrix[x, y] = distance
        distance_matrix[y, x] = distance

# A n x n matrix where n is the number of trajectories. 
pickle.dump(distance_matrix, open("distance_matrix_euclidian.pkl", "wb"))
