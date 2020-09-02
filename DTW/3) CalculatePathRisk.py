import numpy as np
import json
import pickle
from openpyxl import Workbook

# This is a dictionary that maps a given trajectory to two groups. The first is the group of all patients that follow the trajectory
# The second is all patients that follow the trajectory and then go on to experience an opioid event.
file = open('trajectory_groups.pkl', 'rb')
groups = pickle.load(file)

# The groups of clusters is a list of lists. Each sublist contains the index of the trajectories
# in the CUI trajecotry files that are in the cluster
with open('clustering_euclidean_threshold_0.75_two.txt') as outfile:
	clusters = json.load(outfile)

# A dictionarry that mapps a CUI to the initial ICD, CPT, or RxNorm codes that mapped to it.
with open('cui_to_initial2.txt') as outfile:
    mapp = json.load(outfile)

# The trajectores found from the greedy algorithm (of all lengths) after convertion to CUI according to the 
# USMLE metathesaurus.
with open('CUI Trajectories 2.txt') as outfile:
    trajectories = json.load(outfile)

diabetes_075 = [134]

arthropathies_075 = [7,87, 58, 138, 98 ]

normal = 50966/8880896

clusters = sorted(clusters, key = len, reverse=True)
num_clusters = 20

wanted = [0, 1, 4, 6, 7, 8]
wanted2 = [0, 5]


# This function recursively finds all the trajectories of original ICD, CPT, and RxNorm codes that map to the CUI trajectory
def possible_traj(traj):
	trajs = []
	if len(traj) == 0:
		return trajs
	originals = mapp[traj[0]]
	next_level = possible_traj(traj[1:])
	for original in originals:
		if len(next_level) > 0:
			for t in next_level:
				trajs.append([original]+t)
		else:
			trajs.append([original])
	return trajs


counter = 0

# For each cluster of trajectories
for cluster in range(num_clusters):
	# We are trying to find the patients that followed it
	all_patients = set([])
	# And the patients that went on to experience an opioid event.
	opioid = set([])

	# As well as the nodes and edges that make up the trajecotry. 
	nodes = {}
	edge= {}

	cluster_group = clusters[cluster]

	# For the trajectory in the cluster
	for t in cluster_group:

		# Get the trajectory (the trajectories are represented as indecies in the trajectory list)
		traj = trajectories[t]

		# Get the list of possivle trajectories of the CUI trajectory
		converted_trajs = possible_traj(traj)

		# Covert to tuples for dictionary lookup
		converted_trajs = set([tuple(c) for c in converted_trajs])

		# For each of the converted trajectories
		for converted in converted_trajs:
			
			# Get te truncated names of the codes in the trajectory
			converted_names = []
			for node in converted:

				# Either the single letter truncation for ICD
				if len(node) == 3 or node[0:3] == "T40" or node[0:3] == "F11":
					converted_names.append(node[0:1])

				# Opioid for RxNorm
				elif node[0] == "O":
					converted_names.append("opioids")
				# Or other for CPT procedures
				else:
					converted_names.append("other")

			# If the trajectory ends with an opioid event and the trajectory is present in groups (meaning patients follow it)
			if (converted[-1][0:3] == "T40" or converted[-1][0:3] == "F11") and tuple(converted[:-1]) in groups:

				# Add the set of patients that follow it to the running totol of all the patients in the cluster
				all_patients = all_patients | groups[tuple(converted[:-1])]['all']
				opioid = opioid | groups[tuple(converted[:-1])]['opioid']
				
				# Keeps track of the number of patients that pass through each node and edge
				for node in range(len(converted_names)):
					if converted_names[node] not in nodes:
						nodes[converted_names[node]] = groups[tuple(converted[:-1])]['all']
					else:
						nodes[converted_names[node]] = nodes[converted_names[node]] | groups[tuple(converted[:-1])]['all']
					if node > 0:
						if converted_names[node-1] not in edge:
							edge[converted_names[node-1]] = {converted_names[node]: groups[tuple(converted[:-1])]['all']}
						elif converted_names[node] not in edge[converted_names[node-1]]:
							edge[converted_names[node-1]][converted_names[node]] = groups[tuple(converted[:-1])]['all']
						else:
							edge[converted_names[node-1]][converted_names[node]] = edge[converted_names[node-1]][converted_names[node]] | groups[tuple(converted[:-1])]['all']
	counter += 1

	# If the cluster has patients
	if len(all_patients) > 0:

		# And the RR is greater than 1 and at least 10 patients in the cluster have an adverse opioid event
		# Print out all the information fro the cluster
		if (len(opioid)/len(all_patients) / normal) > 1 and len(opioid) >10:
			print(len(all_patients), len(opioid), len(opioid)/len(all_patients) / normal, counter)
			for node in nodes:
				nodes[node] = len(nodes[node])
			print(nodes)
			for node in edge:
				for node2 in edge[node]:
					edge[node][node2] = len(edge[node][node2])
			print(edge)