# Step 6

'''This resturctures the demographics file in a similar way to the other files and removes all patients without
demographic information'''

import json

with open('Demographics to Patient.txt') as outfile:
    data = json.load(outfile)

with open('ICD, Procedures, and Meds Patients to Diagnosis trimmed.txt') as outfile:
	patient_diseases = json.load(outfile)

restructured = {}
for s in data:
	restructured[s] = {}
	for r in data[s]:
		restructured[s][r] = {}
		for a in data[s][r]:
			restructured[s][r][a] = {}
			for patient in data[s][r][a]:
				if patient in patient_diseases:
					restructured[s][r][a][patient] = 0


with open('Demographics to Patient trimmed.txt', 'w') as outfile:
	json.dump(restructured,outfile)

