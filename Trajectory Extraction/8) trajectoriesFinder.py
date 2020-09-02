# Step 8

'''Genrate Trajectories from the ordered pairs'''


import json
import numpy as np

'''A JSON file containing all pairs of ICD code - ICD code relationships and their temporal direction.
{ICD Code1: {ICD Code1a: {"p":p_1a, "RR":RR_1a, "shared":sharedpatients_1b}, ICD Code1b: {"p":p_1b, "RR":RR_1b, "shared":sharedpatients_1b}, ...}, 
ICD Code2: {...}} Each ICD Code maps to a dictionary of all ICD codes that are temporally after it, where each following ICD
code maps to the shared patients between the two'''
with open('ICD, Procedures, and Meds, 10000 All.txt') as outfile:
    data = json.load(outfile)

'''This is a JSON file consisting of a dictionary with the format
{Patient Code1: {ICD Code1: {Date1 : 0, Date2 : 0,...}, {ICD Code2: {Date1: 0}},...},...}
Where each patient maps to all the diagnoses that they received. Each diagnosis also maps to a dictonary
of dates containing all the dates that the patient was diagnosed with that ICD code.
The dates are represented as a integer. Where the number is the number of days since 01/01/0000. This
makes the math easier later and JSONs also don't support datetime'''
with open('ICD, Procedures, and Meds Patients to Diagnosis trimmed.txt') as outfile:
	data2 = json.load(outfile)

'''This is a JSON file consisting of a dictionary with the format
{ICD Code1: {Patient Key1: 0, Patient Key2: 0, ...}, ICD Code2: {Patient Key3: 0,...},...}
Where each ICD code maps to all the patients that recieved that diagnosis at any point
Each Patient Key maps to 0 as I wanted O(1) lookup in the file byt JSON does not support sets.'''
with open('ICD, Procedures, and Meds Diagnosis to Patients trimmed.txt') as outfile:
	data3 = json.load(outfile)

'''Similar to above. Converts the file that takes in a demographic set and returns the patients in that
demographic to return a numpy array.'''
disease_patient = {}
for disease in data3:
	disease_patient[disease] = set(data3[disease].keys())

'''Converts the file of patients to each diagnosis and the discharges of the diagnosis to a format where
the discharge dates are stored in a numpy array'''
patient_disease = {}
for patient in data2:
	patient_disease[patient] = {}
	for diagnosis in data2[patient]:
		patient_disease[patient][diagnosis] = np.array(data2[patient][diagnosis], dtype=float)


paths  = {}

# A function that for a list of diagnoses, finds all the patients that traversed those diagnoses in that order
def shared_patients(diagnoses):
	sets = [disease_patient[diag] for diag in diagnoses]
	patients = set.intersection(*sets)
	total = 0
	for patient in patients:
		dates = np.array([min(patient_disease[patient][diag]) for diag in diagnoses])
		ordered_dates = sorted(dates)
		if numpy.array_equal(ordered_dates, dates):
			total += 1
	return total


trajectories = []
# A recursive function that greedilly finds all paths that at least one patient follows made up of ordered paris. 
def next_level(diagnosis, prior_diagnoses, levels_to_go = 0):
	for d in data[diagnosis]:
		if d not in prior_diagnoses:
			if levels_to_go > 1:
				if d in data:
					next_level(d, prior_diagnoses + [d], levels_to_go - 1)
			else:
				shared_path = shared_patients(prior_diagnoses + [d])
				if shared_path > 0:
					trajectories.append(prior_diagnoses + [d] + [shared_path])

# Calls the function starting with each potiential diagnosis
for d in data:
	# The final number is the desired path length minus one.
	next_level(d, [d], 6)

with open('Trajectories of Length XXXXX.txt', 'w') as outfile:
	json.dump(trajectories,outfile)
