#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import numpy as np
import matplotlib.pyplot as plt
import csv
import random
import collections
import copy
import matplotlib.gridspec as gridspec
from sets import Set
import matplotlib.style
import matplotlib as mpl
mpl.style.use('classic')

SCENARIO='cuboid'

import format as fo

paper_version = True
paper_no_legend = False

# determine if plots are shown
show_plots = True
if len(sys.argv) >= 2:
  show_plots = False
  
remove_outlier = True
outlier_top = 1
outlier_bottom = 2
  
# read csv file
report_filename = "build_release/logs/log.csv"
#report_filename = "build_release/logs/log.0000000.csv"

caption = u'Weak scaling, Hazel Hen'

print("csv file: {}".format(report_filename))
data = []
#data = np.genfromtxt(report_filename, delimiter=';')
#n = np.shape[0]

column_key_map = {}
with open(report_filename) as csvfile:
  reader = csv.reader(csvfile, delimiter=';')
  for i,row in enumerate(reader):
    if len(row) > 0:
      if '#' in row[0]: 
        row[0] = str(row[0][1:]).strip()
        column_keys = row
        for i,key in enumerate(column_keys):
          column_key_map[key] = i
      else:
        data.append(row)
n = len(data)

print("keys: {}".format(column_key_map))

def getCol(colno):
  cols = []
  for i in range(len(data)):
    cols.append(data[i][colno])
  return cols  
  
def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False
    
def isint(value):
  try:
    int(value)
    return True
  except ValueError:
    return False
    
descriptions = {
   0:  "Stamp",
   1:  "Host",
   2:  "NProc",
   3:  "X",
   4:  "Y",
   5:  "Z",
   6:  "F",
   7:  "Total FE",
   8:  "Total M",
   9:  "End Time",
  10:  "Dur. Init",
  11:  "Stretch Sim",
  12:  "Int. Init",
  13:  "Main Sim",
  14:  "Total",
  15:  "Total (User)",
  16:  "Total (System) ",
  17:  "ODE",
  18:  "Parabolic",
  19:  "FE",
  20:  "FE before Main Sim",
  21:  "Mem. Consumption after 1st timestep",
  22:  "Memory Consumption At End ",
  23:  "Parabolic reason",
  24:  "Newton reason",
  25:  "parabolic n. iter",
  26:  "min",
  27:  "max",
  28:  "newton n. iter",
  29:  "min",
  30:  "max",
  31:  "level 0: stimulation handling",
  32:  "level 0: problem solve",
  33:  "level 1: MAIN_TIME_LOOP overhead",
  34:  "level 1: MONODOMAIN_TIME_LOOP overhead",
  35:  "level 1: ELASTICITY_LOOP overhead",
  36:  "level 1: SolverDAE solve",
  37:  "level 1: SolverParabolic solve",
  38:  "level 1: SolverFE solve",
  39:  "level 1: interpolate 1D->3D",
  40:  "level 1: interpolate 3D->1D",
  41:  "level 1: file output",
  42:  "level 2: solver overhead",
  43:  "level 2: 0D solve",
  44:  "level 2: 1D solve",
  45:  "level 2: 3D solve",
  46:  "level 3: 1D assembly",
  47:  "level 3: 1D solve",
  48:  "level 3: 1D other",
  49:  "level 3: 3D assembly",
  50:  "level 3: 3D solve",
  51:  "level 3: 3D other",      
  90:  "FE solver (pre load)",
  91:  "ODE solver (pre load)",
  92:  "parabolic solver (pre load)",
  93:  "file output (pre load)",
  94:  "export EMG (user)",
  95:  "export EMG (system)",
  96:  "file output (user)",
  97:  "file output (system)",
  98:  "file output (system, pre load)",
  99:  "MsolverId",
 100:  "MprecondId",
 101:  "odeSolverId",
 102:  "NumberOfElementsInAtomX",
 103:  "NumberOfElementsInAtomY",
 104:  "NumberOfElementsInAtomZ",
 105:  "nSubdomainsX",
 106:  "nSubdomainsY",
 107:  "nSubdomainsZ",
 108:  "ModelType",
 }

max_index = len(column_key_map)

string_indices = [
  column_key_map["timestamp"],
  column_key_map["hostname"],
  column_key_map["scenarioName"]
]
int_indices = [
  column_key_map["memoryData"],
  column_key_map["memoryPage"],
  column_key_map["memoryResidentSet"],
  column_key_map["memoryVirtual"],
  column_key_map["nInstancesComputedGlobally"],
  column_key_map["nRanks"],
  column_key_map["rankNo"],
  column_key_map["nElements1D"],
  column_key_map["nDofs1D"],
  column_key_map["nNodes1D"],
]

# timestamp;hostname;dt_0D;dt_1D;dt_3D;endTime;memoryData;memoryPage;memoryResidentSet;memoryVirtual;nDofs;nElements;nInstancesComputedGlobally;nNodes;nRanks;rankNo;scenarioName;totalUsertime;duration_0D;n;duration_1D;n;duration_total;n;write output;n;

# all other are float_indices
float_indices = list(set(range(max_index)) - set(int_indices) - set(string_indices))

def extract_data(data):

  datasets = dict()
  
  # extract lines
  n = 0
  for dataset in data:
    
    if len(dataset) < 17:
      print "Warning: invalid data set"
      continue
    
    # copy dataset to new_data
    new_data = dataset
    
    # fill dataset to size of max_index
    if len(dataset) < max_index+1:
      new_data += (max_index+1-len(dataset)) * [0]
      
    # extract some values
    for index in int_indices:
      new_data[index] = int(new_data[index])     if isint(new_data[index])    else 0
    for index in float_indices:
      new_data[index] = float(new_data[index])     if isfloat(new_data[index])    else 0.0
      
    # define sorting key, defines one unique data point (x-axis value without the y-axis value)
      
    scenario_name = column_key_map["scenarioName"]
    if new_data[scenario_name] == "weak_scaling":
      s = "w"
    else:
      s = "s"
      
    nRanks = column_key_map["nRanks"]
    nElements = column_key_map["nElements1D"]
    key = "{}{:05}".format(s,new_data[nRanks])
      
    # store extracted values
    if key not in datasets:
      datasets[key] = dict()
      datasets[key]["value"] = []
    
    datasets[key]['value'].append(new_data)
    n += 1
    
  # compute mean value
  for key in datasets:
    
    result = copy.copy(datasets[key]['value'][0])
    variances = copy.copy(datasets[key]['value'][0])
    result += [0]
    
    for i in int_indices + float_indices:
      
      # reduce values
      result[i] = 0
      value_list = []
      for j in range(len(datasets[key]['value'])):
        value = datasets[key]['value'][j][i]
        
        if value != 0:
          value_list.append(value)
      
      # remove outlier
      value_list = sorted(value_list)
      n = len(value_list)
      
      if n > outlier_bottom+outlier_top and remove_outlier:
        value_list = value_list[outlier_bottom:-outlier_top]
        
      # compute mean and standard deviation
      result[i] = np.mean(value_list)
      variances[i] = np.std(value_list)
      number = n
      #print "mean: {}, std: {}".format(result[i], variances[i])
        
    result[max_index+1] = len(datasets[key]["value"])
        
    datasets[key]['value'] = result
    datasets[key]['variance'] = variances
    datasets[key]['number'] = number
    
  datasets = collections.OrderedDict(sorted(datasets.items()))
    
  return datasets
  
datasets = extract_data(data)

###############################################################
# output to console
print ""
print "------------- duration -------------------------------------------"
print "{:10}, {:7}, {:13}, {:9}, {:4}, {:10}, {:10}, {:10}, {:3}, {:10}".\
format("key", "nproc", "#M", "#M/proc", "#F", "solve: 0D", "1D", "total", "n", "memData")
for key in datasets:
  
  nF = int(datasets[key]["value"][column_key_map["nInstancesComputedGlobally"]])
  nM = int(datasets[key]["value"][column_key_map["nElements1D"]]*datasets[key]["value"][column_key_map["nInstancesComputedGlobally"]])
  nproc = int(datasets[key]["value"][column_key_map["nRanks"]])
  number = datasets[key]["number"]
  
  print("{:10}, {:7}, {:13}, {:9}, {:4}, {:10}, {:10}, {:10}, {:3}, {:10}".\
  format(key, 
  nproc,
  nM,
  nM/nproc,
  nF,
  fo.str_format_seconds(datasets[key]["value"][column_key_map["duration_0D"]]),
  fo.str_format_seconds(datasets[key]["value"][column_key_map["duration_1D"]]),
  fo.str_format_seconds(datasets[key]["value"][column_key_map["duration_total"]]),
  number,
  fo.str_format_memory(datasets[key]["value"][column_key_map["memoryData"]])
  ))
print ""
print ""

###############################################################
#######################################################
# plot
# x-axis: n processes
# y-axis: total time
plt.rcParams.update({'font.size': 20})
plt.rcParams['lines.linewidth'] = 3
plt.rcParams['lines.markersize'] = 8

output_path = ""

colors = {
  column_key_map["duration_total"]: "ko-",      # total
  column_key_map["duration_0D"]: "yd-",      # 0D
  column_key_map["duration_1D"]: "rv-",      # 1D
  "38": "gs-",      # 3D
  "39": "bp-",     # 1D->3D
  "40": "c<-",      # 3D->1D
  "41": "bx-",      # file output
  column_key_map["memoryData"]: "mo-",      # memory consumption  
  
  
  str(column_key_map["duration_total"])+"optimal": "k--",      # total
  str(column_key_map["duration_0D"])+"optimal": "y--",      # 0D
  str(column_key_map["duration_1D"])+"optimal": "r--",      # 1D
  "38optimal": "g--",      # 3D
  
  str(column_key_map["duration_total"])+"efficiency": "ko-",      # total
  str(column_key_map["duration_0D"])+"efficiency": "yd-",      # 0D
  str(column_key_map["duration_1D"])+"efficiency": "rv-",      # 1D
  "38efficiency": "gs-",      # 3D
}

labels = {
  column_key_map["duration_total"]: "total",      # total
  column_key_map["duration_0D"]: "solver 0D model",      # 0D
  column_key_map["duration_1D"]: "solver 1D model",      # 1D
  "38": "solver 3D model",      # 3D
  "39": u"homogenization, 1D to 3D",     # 1D->3D
  "40": u"interpolation, 3D to 1D",      # 3D->1D
  "41": "file output",      # file output
  column_key_map["memoryData"]: "memory consumption",      # memory consumption
}
#plotkeys = [13, 17, 18, 19, 20]

######################
# create plot cubes
caption = "Strong scaling, Hazel Hen, 1e5 1D elements"
outfile = output_path+SCENARIO+'_strong_scaling.pdf'
if paper_no_legend:
  plt.figure("strong scaling", figsize=(10,12))
else:
  plt.figure("strong scaling", figsize=(14,12))

output_path = ""
plotdata = collections.OrderedDict()
xdata = Set()
xtickslist = []
plotkeys = Set()
plotkeys2 = Set()

# key is the initially defined sorting key
for key in datasets:
  
  if "s" not in key:
    continue
  dataset = datasets[key]['value']
  variances = datasets[key]['variance']
  nproc = dataset[column_key_map["nRanks"]]
  nM = dataset[column_key_map["nElements1D"]]*dataset[column_key_map["nInstancesComputedGlobally"]]
    
  xtickslist.append((nM,nproc))
  
  # loop over different curves (e.g. different measurements)
  for plotkey_number in [column_key_map["duration_0D"],column_key_map["duration_1D"],column_key_map["duration_total"]]:
    
    plotkey = plotkey_number
    
    # define x value and y value
    xvalue = nproc
    yvalue = dataset[plotkey_number]
    yvalue_variance = variances[plotkey_number]
      
    if plotkey not in plotdata:
      plotdata[plotkey] = dict()
      plotdata[plotkey]['value'] = collections.OrderedDict()
      plotdata[plotkey]['variance'] = collections.OrderedDict()
      
    plotdata[plotkey]['value'][xvalue] = yvalue
    plotdata[plotkey]['variance'][xvalue] = yvalue_variance
    xdata.add(xvalue)
    plotkeys.add(plotkey)
    
  # loop over different curves (e.g. different measurements)
  for plotkey_number in [column_key_map["duration_0D"],column_key_map["duration_1D"],column_key_map["duration_total"]]:
    
    plotkey = str(plotkey_number)+"optimal"
    
    # define x value and y value
    xvalue = nproc
    
    # optimal value for duration
    optimal = datasets["s00001"]['value'][plotkey_number] / xvalue
    yvalue = optimal
    yvalue_variance = 0
      
    if plotkey not in plotdata:
      plotdata[plotkey] = dict()
      plotdata[plotkey]['value'] = collections.OrderedDict()
      plotdata[plotkey]['variance'] = collections.OrderedDict()
      
    plotdata[plotkey]['value'][xvalue] = yvalue
    plotdata[plotkey]['variance'][xvalue] = yvalue_variance
    xdata.add(xvalue)
    plotkeys.add(plotkey)

    # parallel efficiency
    plotkey = str(plotkey_number)+"efficiency"
    yvalue = optimal / dataset[plotkey_number]
    yvalue_variance = 0
    
    if plotkey not in plotdata:
      plotdata[plotkey] = dict()
      plotdata[plotkey]['value'] = collections.OrderedDict()
      plotdata[plotkey]['variance'] = collections.OrderedDict()
      
    plotdata[plotkey]['value'][xvalue] = yvalue
    plotdata[plotkey]['variance'][xvalue] = yvalue_variance
    plotkeys2.add(plotkey)


gs = gridspec.GridSpec(2,1,height_ratios=[3,2])
plt.subplot(gs[0])

plt.ylabel('Runtime (s)')
#plt.legend(loc='best')
plt.grid(which='both')

# loop over curves and plot data with given label and color
plotkeys = sorted(plotkeys)
for plotkey in plotkeys:
    
  xlist = sorted(plotdata[plotkey]["value"])
  ylist = [item[1] for item in sorted(plotdata[plotkey]["value"].items())]
  yerr = [item[1] for item in sorted(plotdata[plotkey]["variance"].items())]

  label = None
  if plotkey in labels:
    label = labels[plotkey]
  color = ""
  if plotkey in colors:
    color = colors[plotkey]
  linewidth = 3
  if "optimal" in str(plotkey):
    linewidth=1
  plt.errorbar(xlist, ylist, fmt=color, yerr=yerr, label=label, lw=linewidth)
  
ax = plt.gca()
ax.set_xscale('log', basex=10)
ax.set_yscale('log', basey=10) 

plotkeys = []

if paper_no_legend:
  plt.subplots_adjust(bottom=0.12)
else:
  plt.subplots_adjust(right=0.57, top=0.84, bottom=0.12)
  plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., frameon=False)

#ax.set_xticks(np.linspace(000,60000,5))

plt.subplot(gs[1])
ax = plt.gca()
ax.set_xscale('log', basex=10)
#plt.subplot(212)

# loop over curves and plot data with given label and color
plotkeys2 = sorted(plotkeys2)
for plotkey in plotkeys2:
    
  xlist = sorted(plotdata[plotkey]["value"])
  ylist = [item[1] for item in sorted(plotdata[plotkey]["value"].items())]
  yerr = [item[1] for item in sorted(plotdata[plotkey]["variance"].items())]

  label = None
  if plotkey in labels:
    label = labels[plotkey]
  color = ""
  if plotkey in colors:
    color = colors[plotkey]
  plt.errorbar(xlist, ylist, fmt=color, yerr=yerr, label=label)
  
plt.xlabel('Number of processes')
plt.ylabel('Parallel efficiency (-)')
ax = plt.gca()
ax.set_ylim([0,1.0])
#plt.legend(loc='best')
plt.grid(which='both')

#plt.gcf().subplots_adjust(right=0.89)
if not paper_version:
  plt.title(caption, y=1.1)
  plt.tight_layout()
  
#plt.tight_layout()
plt.savefig(outfile)

if show_plots:
  plt.show()
#quit()

