#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy as np
import subprocess
import datetime
import time
from cycler import cycler
import os

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import itertools

if(len(sys.argv)) >= 2:
  input_filename = sys.argv[1]
else:
  input_filename = "logs/log.csv"
#print("filename: {}".format(input_filename))
	
show_plots = True

n_columns = 60

plt.rcParams.update({'font.size': 16})
plt.rcParams['lines.linewidth'] = 2

def load_df(input_filename):

  # determine number of columns
  with open(input_filename) as f:
    line = f.readline()
    if "~nDofsFiber0" in line:
      pos = line.find("~nDofsFiber0")
      column_names = line[0:pos].split(";")
      n_columns = len(column_names)
    else:
      column_names = line.split(";")
      n_columns = len(column_names)

  # load data frame
  df = pd.read_csv(input_filename, sep=';',  warn_bad_lines=True, comment="#", names=column_names, usecols=range(n_columns), engine='python')

  return df

df = load_df(input_filename)

# Info about the data structure
#print("df info:")
#df.info()
#print(df.head())


#print("n rows:")
#print(len(df.index))

# set options for console display
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

def merge_dicts(x, y):
  z = x.copy()   # start with x's keys and values
  z.update(y)    # modifies z with y's keys and values & returns None
  return z

columns_to_plot = ["durationSolve_implicitSolver"]

df_all = df

df_cg_all = df_all[df_all["implicitSolver_solverType"] == "cg"]
df_cg = df_cg_all[df_cg_all["implicitSolver_preconditionerType"] == "none"]
df_cg_boomeramg = df_cg_all[df_cg_all["implicitSolver_preconditionerType"] == "pchypre"]

df_preonly= df_all[df_all["implicitSolver_solverType"] == "preonly"]
df_gamg = df_preonly[df_preonly["implicitSolver_preconditionerType"] == "gamg"]

fig = plt.figure()
ax = plt.axes()

df_plot=pd.DataFrame()

marker = itertools.cycle(('d','^','o', 's')) 
color=itertools.cycle(('r','g','b','k'))
linestyle_list=['--','-.',':','-']
legend_list=[]

solver_cnt = 0

for solver in ["lu", "cg","gmres", "gamg"]:
  if solver == "cg":
    df_solver=df_cg
  else:
    if solver == "gamg":
      df_solver = df_gamg
    else:
      if solver == "cg_boomeramg":
        df_solver = df_cg_boomeramg
      else:
        df_solver = df_all[df_all["implicitSolver_solverType"] == solver]  
        
  means = df_solver.groupby(['~nElements']).mean()  
  errors = df_solver.groupby(['~nElements']).std()  
  
  means.plot(ax=ax,y="nIterationsTotal_implicitSolver",logx=True, logy=True, yerr=errors, linestyle=linestyle_list[solver_cnt],color=next(color), marker=next(marker)) 
  #means.plot(ax=ax,y="durationSolve_implicitSolver",logx=True, logy=True, yerr=errors, linestyle=linestyle_list[solver_cnt],color=next(color), marker=next(marker))
  
  legend_list.append(solver)  
  solver_cnt +=1
        
ax.grid(which='minor')
ax.set_xlabel("Elements")
#ax.set_ylabel("Duration (s)")
ax.set_ylabel("Iterations (s)")
ax.legend(legend_list)

plt.show()
