# Step 7

'''The main portion of analysis. This performs the cohort match for each patient in each ordered pair.'''

import json
import random
import numpy as np
import math
from scipy import sparse
import multiprocessing as mp

np.seterr(divide='ignore', invalid='ignore')

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
    disease_patient = json.load(outfile)

'''A JSON file containing all pairs of ICD code - ICD code relationships and their temporal direction.
{ICD Code1: {ICD Code1a: shared patients 1a, ICD Code1bL shared patients 1b, ...}, ICD Code2: {...}}
Each ICD Code maps to a dictionary of all ICD codes that are temporally after it, where each following ICD
code maps to the shared patients between the two'''
with open('ICD, Procedures, Med Forward Associations trimmed.txt') as outfile:
    data_forward = json.load(outfile)

'''Maps a paitent identifier number to a dictionary of demographic data for the patient.
{Patient Number 1:{'sex': sex, 'gender': gender, 'age':age}, Patient Number 2: {'sex':sex, ...}...}'''
with open('Patient to Demographics.txt') as outfile:
    patient_to_demographics = json.load(outfile)

'''Maps a demographic info set to the patients that are in that demographic
{sex: {gender: {age: {patient1: 0, patient2: 0}, {age2: {...},...} } gender2: {age:...} } sex2: {gender:{age:{...}}}}'''
with open('Demographics to Patient trimmed.txt') as outfile:
    data = json.load(outfile)

'''Converts the file of patients to each diagnosis and the discharges of the diagnosis to a format where
the discharge dates are stored in a numpy array'''
patient_disease = {}
for patient in data2:
	patient_disease[patient] = {}
	for diagnosis in data2[patient]:
		patient_disease[patient][diagnosis] = np.array(data2[patient][diagnosis], dtype=float)

'''Similar to above. Converts the file that takes in a demographic set and returns the patients in that
demographic to return a numpy array.'''
demographics_to_patient = {}
for sex in data:
	if sex not in demographics_to_patient:
		demographics_to_patient[sex] = {}
	for race in data[sex]:
		if race not in demographics_to_patient[sex]:
			demographics_to_patient[sex][race] = {}
		for age in data[sex][race]:
			demographics_to_patient[sex][race][age] = np.array(list((data[sex][race][age].keys())))



# How many trial sets are run
trial_patients = 10000

#The max number of discharges you can have for a given diagnosis
num_discharges = 10000

# The maximum time after the initial diagnosis that the second diagnosis can come for it to be considered significant
# In this case, 5 years, in seconds. 
time_after_initial = (365 * 5 * 24 * 60 * 60)

#Checks if a patient has a given diagnosis
def get_discharges(patient, d2):
	if d2 in patient_disease[patient]:
		return patient
	else:
		return 0

# Vectorizes the function for application to a matrix of patients. 
get_discharges_vec = np.vectorize(get_discharges, otypes = [int])

# Keeps track of each significant association found
sig_associations = {}

# For each diagnosis that has other diagnoses that come after it
for diagnosis in data_forward:

	# If the diagnosis has patients that have demographic information about them
	if diagnosis in disease_patient:

		# Get all the patients with the diagnosis
		patients = disease_patient[diagnosis]

		# Build list of all discharges for the diagnosis
		d1_set = []
		for patient in patients:
			sex, race, age = patient_to_demographics[patient]['sex'], patient_to_demographics[patient]['race'], patient_to_demographics[patient]['age']
			for date in patient_disease[patient][diagnosis]:

				# Makes sure the date wasn't a placeholder date
				if date != 0:
					d1_set.append((sex, race, str(age), float(date), patient))

		#Trims down the number of discharges to the maximum number allowed.
		# 10,000 discharges are used maximum for a given pair. If there are more, a random selection from the total is used. 
		if len(d1_set) > num_discharges:
			d1_set = random.choices(d1_set, k= num_discharges)

		# Build an np.array of all the dates of each discharge
		dates_of_d1 = np.array([[[discharge[3]]] for discharge in d1_set])
		# Build an np.array of all the patient of each discharge
		patients_of_d1 = np.array([[discharge[4]] for discharge in d1_set])

		# Builds an np.array of (x,y) size where x is the number of discharges for the first diagnosis (max 10,000), y is the 
		# number of trial sets that you want (also 10,000). So for each discharge you now have a list of y patients that mach the patients demographics.
		trial_set = np.array([random.choices(demographics_to_patient[discharge[0]][discharge[1]][discharge[2]], k = trial_patients) for discharge in d1_set])

		# For each diagnosis2 that comes after the first diagnosis
		for diagnosis2 in data_forward[diagnosis]:

			# Find the number of patients who were discharged with diagnosis1 that had diagnosis2 at any point
			patient_set_d2 = get_discharges_vec(patients_of_d1, diagnosis2)
			# Uses sparse matrix to save time.
			patient_set_d2 = sparse.coo_matrix(patient_set_d2)
			# The eventual final matrix for the original group of patients that keeps track of whether each patient had diagnosis2
			# in the required timeframe after diagnosis1.
			patient_set_d2_filler = sparse.lil_matrix(np.zeros(patient_set_d2.get_shape()))
			# Iterates through all the patients that had diagnosis2 at some point
			for i,j,patient in zip(patient_set_d2.row, patient_set_d2.col, patient_set_d2.data):
				patient = str(patient)
				# Gets the difference between the date of each discharge of diagnosis2 and the discharge for diagnosis1.
				dates = patient_disease[patient][diagnosis2] - dates_of_d1[i]
				# Finds the discharge of diagnosis2 that came soonest after diagnosis1
				best_date = min(dates[dates >=0], default = time_after_initial * 10)
				# If is is less than the cutoff (5 years) that patient is marked positive.
				patient_set_d2_filler[i,j] = 1 if best_date < time_after_initial else 0

			# The number of patients who were discharged with diagnosis1 who were discharges with diagnosis2 within 5 years. 
			patient_set_d2 = np.sum(patient_set_d2_filler, axis = 1)

			# Find the number of patients who were discharged with diagnosis1 that had diagnosis2 at any point from the trial set.
			trial_set_d2 = get_discharges_vec(trial_set, diagnosis2)
			# Uses sparse matrix to save time.
			trial_set_d2 = sparse.coo_matrix(trial_set_d2)
			# The eventual final matrix for the trial groups of patients that keeps track of whether each patient had diagnosis2
			# in the required timeframe after diagnosis1.
			trial_set_d2_filler = sparse.lil_matrix(np.zeros(trial_set_d2.get_shape()))
			for i,j,patient in zip(trial_set_d2.row, trial_set_d2.col, trial_set_d2.data):
				patient = str(patient)
				# Gets the difference between the date of each discharge of diagnosis2 and the discharge for diagnosis1.
				dates = patient_disease[patient][diagnosis2] - dates_of_d1[i]
				# Finds the discharge of diagnosis2 that came soonest after diagnosis1
				best_date = min(dates[dates >=0], default = time_after_initial*10)
				# If is is less than the cutoff (5 years) that patient is marked positive.
				trial_set_d2_filler[i,j] = 1 if best_date < time_after_initial else 0

			# This is a list of the number of patients in each trial that had a discharge of diagnosis2 within 5 years of the
			# discharge of diagnosis1 for the patient that they were matched to demographically. 
			trial_set_d2 = np.sum(trial_set_d2_filler, axis = 1)

			# Calculate the related rate and the p-value according to the paper
			RR = np.sum(patient_set_d2)/np.average(trial_set_d2)
			p = np.average(np.where(trial_set_d2 > np.sum(patient_set_d2), 1, 0))
			
			# If it is significant keep track of it.
			if p < 0.05 and RR > 1 and not(np.isnan(RR)):
				if diagnosis not in sig_associations:
					sig_associations[diagnosis] = {}
				sig_associations[diagnosis][diagnosis2] = {"p":p, "RR":RR, "shared":data_forward[diagnosis][diagnosis2]}

# Similar format to ICD, Procedures, Med Forward Associations trimmed.txt except each pair contains info on the p-value
# RR, and the number of patients that are in the pair. 
with open('ICD, Procedures, and Meds 10000 ALL.txt', 'w') as outfile:
	json.dump(sig_associations,outfile)