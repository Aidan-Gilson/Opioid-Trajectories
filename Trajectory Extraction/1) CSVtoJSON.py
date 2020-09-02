# STEP 1
'''Converts all CSVs with procedural, medication and diagnostic data to JSON format for easier use later on.
The format of final files is described at the bottom'''


import csv
import json
from datetime import datetime
import time
import requests
import numpy as np


patient_diagnoses = {}
diagnoses_patient = {}

# Files for converting from CPT numbers to plaintext
with open('CPT Numbers.txt') as outfile:
    one = json.load(outfile)

with open('CPT Numbers to Name.txt') as outfile:
    two = json.load(outfile)

cpt = [int(x) for x in one]
value_to_name = {}
for d in two:
	value_to_name[int(d)] = two[d]

time_index = 2
patient_index = 1
diagnosis_index = 6
cpt = np.array(cpt)


# Extracts CPT procedures from CSV to convert to json format
with open('procedure.csv') as csvfile:
	readCSV = csv.reader(csvfile, delimiter=',')
	line = 0
	for diagnoses_event in readCSV:
		if line == 0:
			pass
		else:
			patient = diagnoses_event[patient_index]
			if len(diagnoses_event[time_index]) > 0 and len(diagnoses_event) >= 6:
				diagnoses = diagnoses_event[diagnosis_index]
				if diagnoses_event[time_index][-1] == 'M':
					time_actual = datetime.strptime(diagnoses_event[time_index], '%m/%d/%Y %I:%M:%S %p').timestamp()
				elif diagnoses_event[time_index][-1] == str(0):
					time_actual = datetime.strptime(diagnoses_event[time_index], '%Y-%m-%d %H:%M:%S.%f').timestamp()
				elif diagnoses_event[time_index] is None:
					time_actual = 0
				else:
					time_actual = diagnoses_event[time_index]
				if len(diagnoses) < 7:
					diagnoses = value_to_name[cpt[cpt < int(diagnoses[0:5])].max()]
					if patient in patient_diagnoses:
						if diagnoses in patient_diagnoses[patient]:
							patient_diagnoses[patient][diagnoses][time_actual] = 0
						else:
							patient_diagnoses[patient][diagnoses] = {time_actual: 0}
					else:
						patient_diagnoses[patient] = {diagnoses:{time_actual: 0}}

					if diagnoses not in diagnoses_patient:
						diagnoses_patient[diagnoses] = {}
					diagnoses_patient[diagnoses][patient] = 0
		line += 1

patient_index = 2
diagnoses_index = 4
time_index = 5

# Extracts all ICD-10 information from CSV
with open('diag_full.csv') as csvfile:
	readCSV = csv.reader(csvfile, delimiter=',')
	line_count = 0
	for diagnoses_event in readCSV:
		if line_count == 0:
			pass
		else:
			patient = diagnoses_event[patient_index]
			diagnoses = diagnoses_event[diagnoses_index][0:3]
			if diagnoses == 'T40' or diagnoses == 'F11':
				diagnoses = diagnoses_event[diagnoses_index]
			if isinstance(diagnoses_event[time_index], datetime):
				time_actual = diagnoses_event[time_index].timestamp()
			elif diagnoses_event[time_index] is None or diagnoses_event[time_index] == '':
				time_actual = 0
			else:
				time_actual = datetime.strptime(diagnoses_event[time_index], '%Y-%m-%d').timestamp()
			if patient in patient_diagnoses:
				if diagnoses in patient_diagnoses[patient]:
					patient_diagnoses[patient][diagnoses][time_actual] = 0
				else:
					patient_diagnoses[patient][diagnoses] = {time_actual: 0}
			else:
				patient_diagnoses[patient] = {diagnoses:{time_actual: 0}}

			if diagnoses not in diagnoses_patient:
				diagnoses_patient[diagnoses] = {}
			diagnoses_patient[diagnoses][patient] = 0
		line_count +=1

patient_index = 1
diagnoses_index = 3
time_index = 2

# Extracts all Medication data
with open('opioid_meds.csv') as csvfile:
	readCSV = csv.reader(csvfile, delimiter=',')
	line_count = 0
	for diagnoses_event in readCSV:
		if line_count == 0:
			pass
		else:
			patient = diagnoses_event[patient_index]
			diagnoses = diagnoses_event[diagnoses_index]
			if diagnoses_event[time_index][-1] == 'M':
				time_actual = datetime.strptime(diagnoses_event[time_index], '%m/%d/%Y %I:%M:%S %p').timestamp()
			elif diagnoses_event[time_index][-1] == str(0):
				time_actual = datetime.strptime(diagnoses_event[time_index], '%Y-%m-%d %H:%M:%S.%f').timestamp()
			elif diagnoses_event[time_index] is None:
				time_actual = 0
			else:
				time_actual = diagnoses_event[time_index]
			if patient in patient_diagnoses:
				if diagnoses in patient_diagnoses[patient]:
					patient_diagnoses[patient][diagnoses][time_actual] = 0
				else:
					patient_diagnoses[patient][diagnoses] = {time_actual: 0}
			else:
				patient_diagnoses[patient] = {diagnoses:{time_actual: 0}}

			if diagnoses not in diagnoses_patient:
				diagnoses_patient[diagnoses] = {}
			diagnoses_patient[diagnoses][patient] = 0
		line_count +=1


'''This is a JSON file consisting of a dictionary with the format
{Patient Code1: {ICD Code1: {Date1 : 0, Date2 : 0,...}, {ICD Code2: {Date1: 0}},...},...}
Where each patient maps to all the diagnoses that they received. Each diagnosis also maps to a dictonary
of dates containing all the dates that the patient was diagnosed with that ICD code.
The dates are represented as a integer. Where the number is the number of days since 01/01/0000. This
makes the math easier later and JSONs also don't support datetime'''
with open('ICD, Procedures, and Meds Patients to Diagnosis.txt', 'w') as outfile:
	json.dump(patient_diagnoses,outfile)

'''This is a JSON file consisting of a dictionary with the format
{ICD Code1: {Patient Key1: 0, Patient Key2: 0, ...}, ICD Code2: {Patient Key3: 0,...},...}
Where each ICD code maps to all the patients that recieved that diagnosis at any point
Each Patient Key maps to 0 as I wanted O(1) lookup in the file byt JSON does not support sets.'''
with open('ICD, Procedures, and Meds Diagnosis to Patients.txt', 'w') as outfile:
	json.dump(diagnoses_patient,outfile)