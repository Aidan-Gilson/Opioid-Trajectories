# Step 3

'''Performs the Fisher Exact test to reduce the number of pair counts down before the cohort matching. Done purly to 
speed up the process. Final associations and directions of pairs is found from the cohort match'''

import json
import time
from fractions import Fraction
import scipy.stats as stats

'''This is a JSON file consisting of a dictionary with the format
{ICD Code1: {Patient Key1: 0, Patient Key2: 0, ...}, ICD Code2: {Patient Key3: 0,...},...}
Where each ICD code maps to all the patients that recieved that diagnosis at any point
Each Patient Key maps to 0 as I wanted O(1) lookup in the file byt JSON does not support sets.'''
with open('ICD, Procedures, and Meds Diagnosis to Patients.txt') as outfile:
    data = json.load(outfile)

'''This is a JSON file consisting of a dictionary with the format
{Patient Code1: {ICD Code1: {Date1 : 0, Date2 : 0,...}, {ICD Code2: {Date1: 0}},...},...}
Where each patient maps to all the diagnoses that they received. Each diagnosis also maps to a dictonary
of dates containing all the dates that the patient was diagnosed with that ICD code.
Each Patient Key maps to 0 as I wanted O(1) lookup in the file byt JSON does not support sets.'''
with open('ICD, Procedures, and Meds Patients to Diagnosis.txt') as outfile:
	patient_diseases = json.load(outfile)



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

# A function that performs the Fisher Exact test.
def fish(a,b,c,d):
	one = choose(a+b, a)
	two = choose(c+d, c)
	three = choose(a+b+c+d, a+c)
	return one*two/three

# The number of test is all the pairwise combinations of the ICD, CPT, or RxNorm codes
num_tests = choose(len(data),2)

# Gets a list and number of patients in the trial
patients = set(patient_diseases.keys())
num_patients = len(patients)

# Converts the dictionary of ICD codes to affected patients to use actual sets rather
# than dictionaries mapping to 0
code_patients = {}
for diagnoses in data:
	code_patients[diagnoses] = set(data[diagnoses].keys())


# A list of all the ICD codes that have been used in a diagnosis and the number of them
all_diagnoses = list(code_patients.keys())
num_diagnoses = len(all_diagnoses)
print(num_diagnoses)


# The dictionary that will eventually store all the found significant associations
significant_associations = {}

current = time.time()

#Iterates through the list of diagnoses 
for diagnosis1_index in range(num_diagnoses-1, 0, -1):
	diagnosis1 = all_diagnoses[diagnosis1_index]
	significant_associations[diagnosis1] = {}

	# Keeping track of time
	print(num_diagnoses- diagnosis1_index, time.time()-current, len(code_patients[diagnosis1]))
	current = time.time()

	# Only iterates through diagnoses that have not already been compared to diagnosis1 so
	# test are not performed twice.
	for diagnosis2_index in range(diagnosis1_index):
		diagnosis2 = all_diagnoses[diagnosis2_index]

		# Shouln't be used, but just makes sure diagnosis1 and diagnosis2 are not the same
		if diagnosis1 != diagnosis2:

			# Finds the shared patients between the two diagnoses
			shared_patients = code_patients[diagnosis1] & code_patients[diagnosis2]
			num_shared = len(shared_patients)

			# An arebitrarilly chosen cutoff but I just wanted to make sure there was a 
			# non-trivial amount of shared patients
			if num_shared > 10 or diagnosis1[0:3] == 'T40' or diagnosis1[0:3] == 'F11' or diagnosis2[0:3] == 'T40' or diagnosis2[0:3] == 'F11':

				# Finds the number of patients that have diagnosis1 but not diagnosis2
				d1_positive = len(code_patients[diagnosis1] - shared_patients)

				# Finds the number of patients that have diagnosis2 but not diagnosis1
				d2_positive = len(code_patients[diagnosis2] - shared_patients)

				# Finds the number of patients that have neither diagnosis
				unaffected = num_patients - num_shared

				# Performs the fisher exact test on the two diagnoses
				# p = fish(d1_positive,unaffected,num_shared,d2_positive)
				holder , p = stats.fisher_exact([[num_shared, d1_positive], [d2_positive, unaffected]])

				# If the p-value is less than 0.05 accounting for the Bonferonni correction then...
				if p < 0.05/num_tests:

					# ...saves the associations allong with the p value.
					significant_associations[diagnosis1][diagnosis2] = p

# Dumps the data into a new JSON file.
'''The significant associations found from the FisherExact.py file. Data is in the format:
{Diagnosis1: {Diagnosis1A: p1A, Diagnosis1B: p1B,...}, Diagnosis2: {Diagnosis2A: p2A,...},...}'''
with open('ICD, Procedures, and Meds Fisher Exact.txt', 'w') as outfile:
	json.dump(significant_associations,outfile)


			