#postprocess
import sys
import numpy as np
import csv
import matplotlib.pyplot as plt
from collections import defaultdict
import math
SCENARIO = 'hodgkin-huxley'

log_filename = "build_debug/logs/log.csv"

caption = u'scaling, multigrid solvers'

data = []

column_key_map = {}
with open (log_filename) as csvfile:
	reader=csv.reader (csvfile, delimiter= ';')
	for i,row in enumerate(reader):
		if len(row) > 0:
			if '#' in row[0]:
				row[0] = str(row[0][1:]).strip()
				column_keys = row
				for i,key in enumerate(column_keys):
					column_key_map[key] = id
			else:
				data.append(row)
n = len (data)

## timestamp;hostname;anonymous0_preconditionerType;anonymous0_solverType;anonymous1_preconditionerType;anonymous1_solverType;memoryData;memoryPage;memoryResidentSet;memoryVirtual;nDofs;nElements;nNodes;nRanks;rankNo;scenarioName;totalUsertime;write output;n;
def extract_data(data):
	scenario_data_map = []
	
	for dataset in data:
		index = 0
		new_data = {}
		for entry in dataset:
			if index == 6:
				new_data['memoryData'] = entry
			elif index == 7:
				new_data['memoryPage']= entry
			elif index == 8:
				new_data['memoryResidentSet'] = entry
			elif index == 9:
				new_data['memoryResidentSet'] = entry
			elif index == 10:
				new_data['memoryVirtual'] = entry
			elif index == 11:
				new_data['nDofs'] = entry
			elif index == 12:
				new_data['nElements'] = entry
			elif index == 15:
				new_data['nRanks'] = entry
			elif index == 16:
				new_data['rankNo'] = entry
			elif index == 17:
				new_data['scenarioName'] = entry
			elif index == 18:
				new_data['totalUsertime'] = entry
			elif index == 19:
				new_data['duration_0D'] = entry
			elif index == 21:
				new_data['duration_1D'] = entry
			elif index == 23:
				new_data['duration_total'] = entry
			
			if index > (len(dataset) -1):
				index = 0
			else:
				index += 1
		scenario_data_map.append(new_data)

	return scenario_data_map
	

raw_data = extract_data(data)

raw_data = sorted(raw_data, key = lambda x: ( x["scenarioName"], x['nRanks']))

plot_data = defaultdict(list)
tempY =[]
name = "none"
ranks = 0
raw_data.append({'scenarioName': "Dummy", 'nRanks': '-1', 'totalUsertime' : 0, 'duration_0D' : 0, 'duration_1D' : 0, 'duration_total' : 0})
for dataset in raw_data:
	newScenario = (dataset['scenarioName'] != name) or (int (dataset['nRanks']) != ranks)
	if newScenario:
		if ranks != 0 or int (dataset['nRanks']) == -1:
			plot_data[name].append((sum(tempY) / len (tempY))) 
		ranks = int(dataset['nRanks']) 
		name = dataset['scenarioName']
	tempY.append(float(dataset['totalUsertime']))
###############################################
#create plots
caption = "strong scaling"
plt.figure ("strong scaling", figsize = (8,8))
output_file = SCENARIO + '_strong_scaling.png'
output_path = ""

#weak scaling x  = number of processes
#weak scaling y = seconds

plt.xlabel('nProcesses')
plt.ylabel('computation time seconds')

# 1, 2, 4, 8, 16, 32
xAxis = [1,2,4]

for key in plot_data:
	plt.plot (xAxis, plot_data[key], label=key)

plt.legend()
plt.show()
##################################################
#extract data
#solver : X,Y - values
plot_data = {}
#for dataset in raw_data:
#	dataset['scenarioName'] = (dataset['totalUsertime'], dataset['nRanks']) 
	
