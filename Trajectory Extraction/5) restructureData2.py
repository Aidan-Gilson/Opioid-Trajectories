# Step 5

'''Not necessary for analysis, but done simply to speed up the cohort analysis. This file resturctures the main JSON files
so that they contain only pairs relevent to adverse opioid events. It first removes Z and R blocks in the ICD system.
It also performs a BFS search from all advers opioid events and finds only the codes that can be reached from them by traversing
an ordered pair (forward or backwards just to be safe). This prevents pairs from being checked later that have no chance of
being in a trajectory.'''


import json
import numpy as np
import math
import datetime

# Loads in the previous files
with open('ICD, Procedures, and Meds Patients to Diagnosis.txt') as outfile:
	patient_diseases = json.load(outfile)

with open('ICD, Procedures, and Meds Diagnosis to Patients.txt') as outfile:
    disease_patient = json.load(outfile)

with open('ICD, Procedures, Med Forward Associations.txt') as outfile:
    data_forward = json.load(outfile)

with open('ICD, Procedures, Med Reversed Associations.txt') as outfile:
    data_reversed = json.load(outfile)

with open('Patient to Demographics.txt') as outfile:
    patient_to_demographics = json.load(outfile)

# These are the possible values of each demographic type that can be used later in the analysis. Any patient that does not have
# a value contained in these will be excluded.
sexes = {'Male', 'Female'}

races = {'Asian', 'American Indian or Alaska Native', 'Black', 'White', 'Other', 'Native Hawaiian or Other Pacific Islander','Unknown'}

ages = [10*x for x in range(0,12)]


# ICD codes corresponding to adverse opioid events
relavent_diagnoses = ['T40.2X4A', 'T40.3X5A', 'T40.601S', 'T40.692D', 'F11.14', 'T40.2X5S', 'T40.604A', 'T40.4X2D', 'T40.3X4S', 'T40.692A', 'T40.4X2S', 'F11.23', 'T40.3X1D', 
'T40.605A', 'F11.90', 'F11.122', 'F11.250', 'F11.982', 'T40.4X5A', 'F11.150', 'T40.3X2A', 'T40.604D', 'F11.188', 'T40.603A', 'T40.1X3D', 'F11.221', 'T40.0X5S', 'T40.3X2D', 
'T40.605S', 'T40.693D', 'T40.0X4A', 'T40.604S', 'F11.120', 'F11.229', 'F11.251', 'F11.93', 'T40.4X4D', 'T40.1X2S', 'T40.0X4S', 'F11.129', 'F11.929', 'T40.3X1S', 'T40.2X2A', 
'F11.29', 'T40.3X2S', 'T40.694S', 'F11.281', 'T40.695S', 'T40.0X2D', 'T40.605D', 'T40.2X3D', 'T40.4X2A', 'T40.2X2D', 'T40.0X3S', 'F11.10', 'T40.4X4A', 'T40.2X1D', 'F11.151', 
'T40.3X3S', 'T40.0X1A', 'T40.692S', 'T40.4X3A', 'T40.3X3A', 'T40.4X3S', 'T40.0X2S', 'T40.2X5D', 'F11.950', 'F11.288', 'T40.2X3A', 'T40.602D', 'T40.602S', 'F11.19', 'T40.2X1S',
'T40.0X5A', 'T40.3X4D', 'T40.4X3D', 'T40.693S', 'T40.1X3S', 'F11.99', 'T40.603S', 'T40.1X4D', 'T40.3X3D', 'F11.959', 'F11.24', 'T40.2X3S', 'T40.0X3D', 'T40.691A', 'F11.920', 
'T40.602A', 'T40.2X4D', 'T40.1X4S', 'F11.988', 'T40.1X2D', 'F11.182', 'T40.695A', 'T40.4X5S', 'T40.1X1A', 'F11.94', 'T40.2X1A', 'T40.691D', 'F11.921', 'F11.951', 'F11.181', 
'T40.3X5D', 'F11.121', 'T40.2X5A', 'T40.603D', 'T40.3X4A', 'T40.0X1D', 'T40.0X5D', 'T40.1X2A', 'F11.20', 'T40.1X3A', 'T40.0X3A', 'F11.282', 'F11.222', 'T40.4X4S', 'F11.259', 
'T40.601A', 'T40.3X5S', 'T40.2X2S', 'T40.0X4D', 'T40.4X5D', 'F11.220', 'F11.922', 'T40.3X1A', 'T40.4X1S', 'T40.1X1S', 'T40.2X4S', 'T40.694A', 'F11.981', 'T40.1X1D', 'T40.4X1D',
 'T40.0X2A', 'T40.694D', 'T40.695D', 'T40.693A', 'F11.159', 'T40.1X4A', 'T40.4X1A', 'T40.691S', 'T40.0X1S', 'T40.601D']

# Never used, but this was a cutoff that could be used to remove pairs with less than the cutoff amount of patients between them
# Pairs with very few patients are almost always excluded anyways because the do not reach any level of statistical significance. 
cutoff = 0

# Just a BFS function. 
def bfs(graph, q):

	queue = q.copy()
	original = queue.copy()

	visited = set(queue.copy())
	relavent = set()

	while queue:
		node = queue.pop(0)

		try:
			neighbors = graph[node].keys()
			if node in original and len(neighbors) > 0:
				relavent.add(node)
		except:
			neighbors = []


		for node2 in neighbors:
			if node2 not in visited:
				visited.add(node2)
				queue.append(node2)
				if node2 not in original:
					relavent.add(node2)
	return relavent

# This removes pairs containing a R or Z code.
data_forward_trimmed = {}
for diagnosis in data_forward:
	if (diagnosis[0:1] != 'Z' and diagnosis[0:1] != 'R' and data_forward[diagnosis] != {}):
		for diagnosis2 in data_forward[diagnosis]:
			if (diagnosis2[0:1] != 'Z' and diagnosis2[0:1] != 'R'):
				if diagnosis not in data_forward_trimmed:
					data_forward_trimmed[diagnosis] = {}
				if data_forward[diagnosis][diagnosis2] >= cutoff:
					data_forward_trimmed[diagnosis][diagnosis2] = data_forward[diagnosis][diagnosis2]

data_reversed_trimmed = {}
for diagnosis in data_reversed:
	if (diagnosis[0:1] != 'Z' and diagnosis[0:1] != 'R' and data_reversed[diagnosis] != {}):
		for diagnosis2 in data_reversed[diagnosis]:
			if (diagnosis2[0:1] != 'Z' and diagnosis2[0:1] != 'R'):
				if diagnosis not in data_reversed_trimmed:
					data_reversed_trimmed[diagnosis] = {}
				if data_reversed[diagnosis][diagnosis2] >= cutoff:
					data_reversed_trimmed[diagnosis][diagnosis2] = data_reversed[diagnosis][diagnosis2]

# This finds all the nodes that can either be reached from an adverse opioid event, or that can reach an adverse opioid event
# through some path of ordered pairs. 
relavent = bfs(data_forward_trimmed, relavent_diagnoses)
relavent1 = set(relavent.copy())
relavent = bfs(data_reversed_trimmed, relavent_diagnoses)
relavent = list(relavent)
relavent_total = set(relavent.copy()) or relavent1


# Creates a dictionary like ICD, Procedures, Med Forward Associations.txt off all pairs containing reachable nodes from the BFS.
reachable = {}
for diag in data_forward_trimmed:
	if diag in relavent_total:
		for diag2 in data_forward_trimmed[diag]:
			if diag2 in relavent_total:
				if diag not in reachable:
					reachable[diag] = {}
				reachable[diag][diag2] = data_forward_trimmed[diag][diag2]


for diag in data_reversed_trimmed:
	if diag in relavent_total:
		for diag2 in data_reversed_trimmed[diag]:
			if diag2 in relavent_total:
				if diag2 not in reachable:
					reachable[diag2] = {}
				reachable[diag2][diag] = data_reversed_trimmed[diag][diag2]


# Restructures the main patient to diagnoses and diagnoses to patient files to contain only patients who have proper demographic
# information and only diagnoses that are reachable from an adverse opioid event. 
patient_diseases_restructured1 = {}
for patient in patient_diseases:
	if patient in patient_to_demographics:
		sex, race, age = patient_to_demographics[patient]['sex'], patient_to_demographics[patient]['race'], patient_to_demographics[patient]['age']
		if sex in sexes and race in races and age in ages:
			for disease in patient_diseases[patient]:
				if disease in relavent_total:
					dates_list = []
					dates = patient_diseases[patient][disease]
					for date in dates:
						if date != 0 and date != '':
							dates_list.append(date)
					if len(dates_list) != 0:
						if patient not in patient_diseases_restructured1:
							patient_diseases_restructured1[patient] = {}
						patient_diseases_restructured1[patient][disease] = dates_list

disease_patient_restructured = {}
for disease in disease_patient:
	if disease in relavent_total:
		for patient in disease_patient[disease]:
			if patient in patient_diseases_restructured1:
				if disease in patient_diseases_restructured1[patient]:
					if disease not in disease_patient_restructured:
						disease_patient_restructured[disease] = {}
					disease_patient_restructured[disease][patient] = 0


# These files are identical in structure to their predicessors, just smaller.
with open('ICD, Procedures, and Meds Patients to Diagnosis trimmed.txt', 'w') as outfile:
	json.dump(patient_diseases_restructured1,outfile)

with open('ICD, Procedures, and Meds Diagnosis to Patients trimmed.txt', 'w') as outfile:
	json.dump(disease_patient_restructured,outfile)

with open('ICD, Procedures, Med Forward Associations trimmed.txt', 'w') as outfile:
	json.dump(reachable,outfile)


