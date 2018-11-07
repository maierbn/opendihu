#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import numpy as np
#import matplotlib.pyplot as plt
import csv
import random
import collections
import copy
#from sets import Set
#import matplotlib.style
#import matplotlib as mpl
#mpl.style.use('classic')

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

caption = u'Performance measurements, Hazel Hen'

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
   
max_index = len(column_key_map)

string_indices = [
  column_key_map["timestamp"],
  column_key_map["hostname"],
  column_key_map["scenarioName"],
  column_key_map["implicitSolver_preconditionerType"],
  column_key_map["implicitSolver_solverType"],
]
int_indices = [
  column_key_map["memoryData"],
  column_key_map["memoryPage"],
  column_key_map["memoryResidentSet"],
  column_key_map["memoryVirtual"],
  column_key_map["nInstancesComputedGlobally"],
  column_key_map["nRanks"],
  column_key_map["rankNo"],
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
      print("Warning: invalid data set")
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
      
    key_solver = column_key_map["implicitSolver_solverType"]
    nRanks = column_key_map["nRanks"]
    scenario_name = column_key_map["scenarioName"]
    solver_type = new_data[key_solver]
    
    key = "{:30}".format(new_data[scenario_name])
      
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
print("")
print("------------- duration -------------------------------------------")
print("{:30}, {:7},  {:10}, {:10}, {:10}, {:3}, {:10}".\
format("key", "nproc", "solve: 0D", "1D", "total", "n", "memData"))
for key in datasets:
  
  nF = int(datasets[key]["value"][column_key_map["nInstancesComputedGlobally"]])
  nproc = int(datasets[key]["value"][column_key_map["nRanks"]])
  number = datasets[key]["number"]
  
  print("{:30}, {:7}, {:10}, {:10}, {:10}, {:3}, {:10}".\
  format(key, 
  nproc,
  fo.str_format_seconds(datasets[key]["value"][column_key_map["duration_0D"]]),
  fo.str_format_seconds(datasets[key]["value"][column_key_map["duration_1D"]]),
  fo.str_format_seconds(datasets[key]["value"][column_key_map["duration_total"]]),
  number,
  fo.str_format_memory(datasets[key]["value"][column_key_map["memoryData"]])
  ))
print("")
print("")
sys.exit(0)  
###############################################################
#######################################################
# plot
# x-axis: n processes
# y-axis: total time
plt.rcParams.update({'font.size': 16})
plt.rcParams['lines.linewidth'] = 3
plt.rcParams['lines.markersize'] = 8

output_path = ""

colors = {
  column_key_map["duration_total"]: "ko-",      # total
  column_key_map["duration_0D"]: "yd-",      # 0D
  "cg_"+str(column_key_map["duration_1D"]): "rv--",      # 1D CG
  "lu_"+str(column_key_map["duration_1D"]): "gd-.",      # 1D LU
  "gmres_"+str(column_key_map["duration_1D"]): "bo:",      # 1D GMRES
  "gamg_"+str(column_key_map["duration_1D"]): "cs--",      # 1D GAMG
  "38o": "gs-",      # 3D
  "39o": "bp-",     # 1D->3D
  "40o": "c<-",      # 3D->1D
  "41o": "bx-",      # file output
  column_key_map["memoryData"]: "mo-",      # memory consumption
  "22o": "mo-",      # memory consumption
  "ghosto": "ko-",   # ghost layer size
  
  "15-": "ko--",      # total
  "36-": "yd--",      # 0D
  "37-": "rv--",      # 1D
  "38-": "gs--",      # 3D
  "39-": "bp--",     # 1D->3D
  "40-": "c<--",      # 3D->1D
  "41-": "bx-",      # file output
  "22-": "mo--",      # memory consumption
  
  "ghost-": "ko--",   # ghost layer size
}

labels = {
  "cg_"+str(column_key_map["duration_total"]): "total (CG)",      # total
  "cg_"+str(column_key_map["duration_0D"]): "solver 0D model (CG)",      # 0D
  "cg_"+str(column_key_map["duration_1D"]): "solver 1D model (CG)",      # 1D
  "lu_"+str(column_key_map["duration_1D"]): "solver 1D model (LU)",      # 1D
  "gmres_"+str(column_key_map["duration_1D"]): "solver 1D model (GMRES)",      # 1D
  "gamg_"+str(column_key_map["duration_1D"]): "solver 1D model (AMG)",      # 1D
  "38o": "solver 3D model",      # 3D
  "39o": u"homogenization, 1D to 3D",     # 1D->3D
  "40o": u"interpolation, 3D to 1D",      # 3D->1D
  "41o": "file output",      # file output
  column_key_map["memoryData"]: "memory consumption",      # memory consumption
  "22o": "mem. cubes",      # memory consumption
  "22-": "mem. elongated",      # memory consumption
  "ghosto": "ghost layer size"   # ghost layer size
}
#plotkeys = [13, 17, 18, 19, 20]

######################
# create plot solvers
caption = "Multi-node weak scaling solvers, Hazel Hen,\n100 el./core"
outfile = output_path+SCENARIO+'_solver_scaling.pdf'
if paper_no_legend:
  plt.figure("linear solvers weak scaling", figsize=(10,12))
else:
  plt.figure("linear solvers weak scaling", figsize=(14,8))

output_path = ""
plotdata = collections.OrderedDict()
xdata = Set()
xtickslist = []
plotkeys = Set()

# key is the initially defined sorting key
for key in datasets:
  
  dataset = datasets[key]['value']
  variances = datasets[key]['variance']
  nproc = dataset[column_key_map["nRanks"]]
  nM = dataset[column_key_map["nElements1D"]]*dataset[column_key_map["nInstancesComputedGlobally"]]
  
  xtickslist.append((nM,nproc))
  
  # loop over different curves (e.g. different measurements)
  for plotkey_number in [column_key_map["duration_1D"]]:
    
    plotkey = dataset[column_key_map["implicitSolver_solverType"]]+"_"+str(plotkey_number)
    
    # define x value and y value
    xvalue = nM
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
  print ("label:",label,", color:",color)
  plt.errorbar(xlist, ylist, fmt=color, yerr=yerr, label=label)
  
  
ax = plt.gca()
ax.set_xscale('log', basex=10) 
ax.set_yscale('log', basey=10)
#ax.set_xscale('log', basey=2) 
#ticks = list(np.linspace(10**4, 10**5, 10)) + list(np.linspace(10**5, 10**6, 10))
#ax.set_xticks(ticks)
#ax.set_xticklabels([int(i/1000.) for i in ticks])
#ax.set_xticks(np.linspace(000,60000,5))

plt.xlabel('Number of 1D elements')
plt.ylabel('Runtime (s)')
ax.set_ylim(0, ax.get_ylim()[1])
#plt.legend(loc='best')
plt.grid(which='major')

if not paper_no_legend:
  print ("legend")
  if paper_version:
    plt.subplots_adjust(right=0.58, top=0.84)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
  else:
    plt.legend(loc='best')

# twin axes for processes
ax2 = ax.twiny()
ax2.set_xlim(ax.get_xlim())
ax2.set_xscale('log', basex=10)

xtickslist = sorted(list(set(xtickslist)))
#xtickslist = [(item[0],int(np.ceil(item[1]/24.))) for item in xtickslist]  # number of nodes instead number of processes

# only leave certain values for number of processes
#show_processes = [1, 2, 3, 4, 6, 8, 16, 32]
show_processes = [1, 3, 10, 32, 100, 316, 1000, 3162, 10000]
xtickslist_new = list()
for item in xtickslist:
  
  # omit number of process values that are already present in xtickslist_new
  if item[1] in [itema[1] for itema in xtickslist_new]:
    continue
  
  # if the current number of process value is in show_processes, add it to xtickslist_new
  if item[1] in show_processes:
    xtickslist_new.append(item)

#xtickslist = list(xtickslist_new)
xticks = [item[0] for item in xtickslist]
xlabels = [int(item[1]) for item in xtickslist]
#xlabels = [int(np.ceil(item[1]/24.)) for item in xtickslist]
print (xtickslist)
print ("xticks:",xticks,", xlabels:",xlabels)

ax2.set_xticks(xticks)
ax2.set_xticklabels(xlabels)
ax2.set_xlabel(r"Number of processes")

if not paper_version:
  plt.title(caption, y=1.1)
  plt.tight_layout()
  
print("write file '{}'".format(outfile))
plt.savefig(outfile)

if show_plots:
  plt.show()
#quit()

