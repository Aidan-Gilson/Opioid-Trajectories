# Step 2
'''Extracts the demographic information from the CSV for all the patients who had at least one encounter in the Yale Healthcare system'''


import csv
import json
import math

patient_to_demographics = {}
demographics_to_patient = {}

races = set()
sexes = set()
ages = set()

# Groups the ages into decades of life for cohort matching later
def roundDown(age):
	try:
		return int(math.floor(int(age)//10)) * 10
	except:
		return(age)

with open('ICD, Procedures, and Meds Patients to Diagnosis.txt') as outfile:
	patient_diseases = json.load(outfile)

with open('cohort_redo2.csv') as csvfile:
	readCSV = csv.reader(csvfile, delimiter=',')
	line = 0
	for patient1 in readCSV:
		if line > 0:
			try:
				patient, age, race, sex = patient1[2], patient1[4], patient1[10], patient1[9]
				roundedAge = roundDown(age)
				if race == 'Other/Not Listed':
					race = 'Other'
				if race == 'White or Caucasian':
					race = 'White'
				if race == 'Black or African American':
					race = 'Black'
				if race == 'Other Pacific Islander' or race == 'Native Hawaiian':
					race = 'Native Hawaiian or Other Pacific Islander'
				if race == '' or race == 'Patient Refused' or race == '*Not Applicable':
					race = 'Unknown'
				if age != '':
					if patient in patient_diseases:
						if patient not in patient_to_demographics:
							patient_to_demographics[patient] = {'sex':sex, 'race':race, 'age':roundedAge} 
						if sex not in demographics_to_patient:
							sexes.add(sex)
							demographics_to_patient[sex] = {race:{roundedAge:{patient:0}}}
						elif race not in demographics_to_patient[sex]:
							races.add(race)
							demographics_to_patient[sex][race] = {roundedAge:{patient:0}}
						elif roundedAge not in demographics_to_patient[sex][race]:
							ages.add(roundedAge)
							demographics_to_patient[sex][race][roundedAge] = {patient:0}
						else:
							demographics_to_patient[sex][race][roundedAge][patient] = 0
			except:
				print(patient)
		else:
			line = 1

'''Maps a paitent identifier number to a dictionary of demographic data for the patient.
{Patient Number 1:{'sex': sex, 'gender': gender, 'age':age}, Patient Number 2: {'sex':sex, ...}...}'''
with open('Patient to Demographics.txt', 'w') as outfile:
	json.dump(patient_to_demographics,outfile)

'''Maps a demographic info set to the patients that are in that demographic
{sex: {gender: {age: {patient1: 0, patient2: 0}, {age2: {...},...} } gender2: {age:...} } sex2: {gender:{age:{...}}}}'''
with open('Demographics to Patient.txt', 'w') as outfile:
	json.dump(demographics_to_patient,outfile)