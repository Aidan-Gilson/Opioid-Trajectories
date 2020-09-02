# Step 4

'''Performs the binomial test to reduce the number of pair counts down before the cohort matching. This takes pairs
from the fisher exact test and reduces them further by testing for preliminary directionallity. Done purly to 
speed up the process. Final associations and directions of pairs is found from the cohort match'''

import json
import datetime
from fractions import Fraction
from decimal import Decimal

'''This is a JSON file consisting of a dictionary with the format
{ICD Code1: {Patient Key1: 0, Patient Key2: 0, ...}, ICD Code2: {Patient Key3: 0,...},...}
Where each ICD code maps to all the patients that recieved that diagnosis at any point
Each Patient Key maps to 0 as I wanted O(1) lookup in the file byt JSON does not support sets.'''
with open('ICD, Procedures, and Meds Diagnosis to Patients.txt') as outfile:
    data = json.load(outfile)

# Converts the dictionary of ICD codes to affected patients to use actual sets rather
# than dictionaries mapping to 0
code_patients = {}
for diagnoses in data:
	code_patients[diagnoses] = set(data[diagnoses].keys())

'''This is a JSON file consisting of a dictionary with the format
{Patient Code1: {ICD Code1: {Date1 : 0, Date2 : 0,...}, {ICD Code2: {Date1: 0}},...},...}
Where each patient maps to all the diagnoses that they received. Each diagnosis also maps to a dictonary
of dates containing all the dates that the patient was diagnosed with that ICD code.
Each Patient Key maps to 0 as I wanted O(1) lookup in the file byt JSON does not support sets.'''
with open('ICD, Procedures, and Meds Patients to Diagnosis.txt') as outfile:
	patient_diseases = json.load(outfile)

'''The significant associations found from the FisherExact.py file. Data is in the format:
{Diagnosis1: {Diagnosis1A: p1A, Diagnosis1B: p1B,...}, Diagnosis2: {Diagnosis2A: p2A,...},...}'''
with open('ICD, Procedures, and Meds Fisher Exact.txt') as outfile:
	significant_associations = json.load(outfile)

# Given a patient and two diagnoses this function returns which diagnosis was first diagnosed in the patient
# The function returns 1 if diagnosis1 was diagnosed first, 0 if diagnosis2 was diagnosed first and 2 if 
# they were diagnosed at the same time or it was unable to be determined.
def timeTest(patient, disease1, disease2):

	# Gets all of the dates that a patient was diagnosed with diagnosis1
	# The exception was for cases in which the patient was diagnosed but the date
	# was unclear and so therefore not stored
	try:
		dates1 = patient_diseases[patient][disease1]
	except:
		print(patient, disease1, disease2)
		return 2

	# Gets all of the dates that a patient was diagnosed with diagnosis2
	try:
		dates2 = patient_diseases[patient][disease2]
	except:
		print(patient, disease1, disease2)
		return 2


	dates1_list = []
	dates2_list = []

	# Creates a list of all the dates that a patient had diagnosis1
	for date in dates1:

		# This removes some data were the date was either empty or it was filled with a 
		# placeholder value of 01/01/1900 when the JSON files were made from the CSV
		if date != 0 and date != '':

			# This converts each date to a datetime object for easy comparison. JSON does not support 
			# datetime so they could not be stored as such.
			dates1_list.append(date)

	# This does the same thing for diagnosis2
	for date in dates2:
		if date != 0 and date != '':
			dates2_list.append(date)

	# Sorts each list of dates
	dates1_list.sort()
	dates2_list.sort()

	# Compares the first element in each list to see which is earlier. There were some exceptions
	# where one/both of the lists was empty so that is why it is in a try/except block.
	try:
		# If diagnosis1 was earlier return 1
		if dates1_list[0] < dates2_list[0]:
			return 1
		# If diagnosis2 was earlier return 0
		elif dates1_list[0] > dates2_list[0]:
			return 0 
		# If they were at the same time return 2
		else:
			return 2
	except:
		return 2

# An n choose k funciton
def choose(n,k):
	if k > n//2:
		k = n - k

	# Had to use Fraction and Decimal imports later on to get around overflow
	# errors if I had been using floats.
	p = Fraction(1)
	for i in range(1,k+1):
		p *= Fraction(n - i + 1, i)
	return int(p)

'''This will store the data such that each diagnosis maps to a dictonary whoses keys are 
diagnoses that were shown to come after diagnosis1. The value of each key will be the number of
shared patients between that diagnosis and the diagnosis mapping to the set of keys i.e.
{Diagnosis1: {Diagnosis2: Shared1-2, Diagnosis3: Shared1-3,...}, Diagnosis2: {Diagnosis3: Shared2-3,...},...}'''
directional_associations = {}

'''This is the same structure as above except each diagnosis maps to a set of keys that were shown 
to come before the original diagnosis, insead of after. Each key still maps to the shared patients.'''
reversed_directional_associations = {}

# Used to track the progress of the analysis
total = 0

# Iterates through each of the diagnoses for their significant associations
for disease1 in significant_associations:

	print("total:", total)

	# Iterates through each diagnosis that was shown to be significantly associated with
	# the first diagnosis. This is frequently empty.
	for disease2 in significant_associations[disease1]:

		# Trackers that keep track of the ammount of patient were diagnosis1 came before...
		before = 0

		# ...and after diagnosis2
		after = 0

		# Finds the shared patients of the diagnoses
		shared_patients = code_patients[disease1] & code_patients[disease2]

		# Iterates through all patients
		for patient in shared_patients:

			# Determines which diagnosis came first and iterates before and after accordingly
			relative_dates = timeTest(patient, disease1, disease2)
			if relative_dates == 1:
				before += 1
			elif relative_dates == 0:
				after += 1

		print(len(shared_patients), before, after)

		# Runs a Bernoulli trial to determine the p value of diagnosis1 coming before2 as many times as it did.
		# Used Decimal to get around float overflow errors from the very large/small numbers.
		p = Decimal(choose(len(shared_patients), before)) * Decimal(0.5) ** Decimal(before) * Decimal(0.5) ** Decimal(after)

		# If is was significant
		if p < 0.05:

			# If diagnosis1 came first store the data in the corresponding way.
			if before > after:
				if disease1 not in directional_associations:
					directional_associations[disease1] = {}
				if disease2 not in reversed_directional_associations:
					reversed_directional_associations[disease2] = {}
				directional_associations[disease1][disease2] = len(shared_patients)
				reversed_directional_associations[disease2][disease1] = len(shared_patients)

			# If diagnosis2 came first store the data in the corresponding way.
			elif after > before:
				if disease2 not in directional_associations:
					directional_associations[disease2] = {}
				if disease1 not in reversed_directional_associations:
					reversed_directional_associations[disease1] = {}
				directional_associations[disease2][disease1] = len(shared_patients)
				reversed_directional_associations[disease1][disease2] = len(shared_patients)
	total += 1


# Store both directions of the temporal associations.
'''A JSON file containing all pairs of ICD code - ICD code relationships and their temporal direction.
{ICD Code1: {ICD Code1a: shared patients 1a, ICD Code1b: shared patients 1b, ...}, ICD Code2: {...}}
Each ICD Code maps to a dictionary of all ICD codes that are temporally after it, where each following ICD
code maps to the shared patients between the two'''
with open('ICD, Procedures, Med Forward Associations.txt', 'w') as outfile:
	json.dump(directional_associations,outfile)
with open('ICD, Procedures, Med Reversed Associations.txt', 'w') as outfile:
	json.dump(reversed_directional_associations,outfile)

