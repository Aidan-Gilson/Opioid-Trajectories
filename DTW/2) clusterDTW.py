# DTW Step 2

'''Clusters the trajecotries based on the pairwise distances found before.'''

import json
import numpy as np
import time
import pickle

# The trajectores found from the greedy algorithm (of all lengths) after convertion to CUI according to the 
# USMLE metathesaurus.
with open('CUI Trajectories 2.txt') as outfile:
	trajectories = json.load(outfile)

# A dictionarry that mapps a CUI to the initial ICD, CPT, or RxNorm codes that mapped to it.
with open('cui_to_initial2.txt') as outfile:
	mapp = json.load(outfile)

# A n x n matrix where n is the number of trajectories. 
file = open('distance_matrix_euclidian_two.pkl', 'rb')
distance_matrix = pickle.load(file)

# The threshold value that a trajectory must be within in order to be added to a group
threshold = 0.75

# Will be a list of lists where each list is a cluster of trajectories.
clusters = []

for traj in range(len(trajectories)):

	# If its the first trajectory, just add it to a cluster.
	if len(clusters) == 0:
		clusters.append([traj])
	else:
		# Keeps track of the minimum distance to a cluster that the trajectory has found
		distance = np.inf
		# The index of the corresponding cluster
		index = None

		# Find the average distance of this trajectory to every trajectory in the cluster.
		for cluster in range(len(clusters)):
			total = 0
			for member in clusters[cluster]:
				total += distance_matrix[traj][member]

			# If its less than the best so far, keep track of the distance and the cluster
			if total/len(clusters[cluster]) < distance:
				distance = total/len(clusters[cluster])
				index = cluster

		# If the best distance is less than the threshold, then add the trajectory to that cluster, else create a new cluster.
		if distance <= threshold:
			clusters[index].append(traj)
		else:
			clusters.append([traj])

# The groups of clusters is a list of lists. Each sublist contains the index of the trajectories
# in the CUI trajecotry files that are in the cluster
with open('clustering_euclidean_threshold_0.75_two.txt', 'w') as outfile:
	json.dump(clusters,outfile)