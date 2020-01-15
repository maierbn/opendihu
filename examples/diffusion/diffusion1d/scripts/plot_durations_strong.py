#!/usr/bin/env python
# -*- coding: utf-8 -*-

# parallel weak scaling

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
#if len(sys.argv) > 1:
#  show_plots = False;

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
fig=plt.figure()
ax2=plt.axes()

df_plot=pd.DataFrame()
marker = itertools.cycle(('o', '+', 's','^','d')) 
color=itertools.cycle(('r','r','r','r','r','g','g','g','g','g','b','b','b','b','b'))

#for solver in ["lu", "cg","gmres", "gamg", "cg_boomeramg"]:
#for solver in ["lu", "cg", "gamg"]:
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
  means.plot(ax=ax,y="durationSolve_implicitSolver",logx=True, logy=True, marker=next(marker),yerr=errors)  
  means.plot(ax=ax2,y="nIterationsTotal_implicitSolver",logx=True, logy=True, marker=next(marker)) 
     
  #for elem in [100,200,500,1000,1480]:
    #df_plot=df_solver[df_solver['~nElements']== elem]
    #means = df_plot.groupby(['nRanks']).mean()
    #errors = df_plot.groupby(['nRanks']).std()
    ##print(means)
    #means.plot(ax=ax,y="durationSolve_implicitSolver",logx=False, logy=True, color=next(color), marker=next(marker), yerr=errors)


#  if df_plot.empty:
#    df_plot=means
#  else:
#    df_plot=pd.DataFrame.append(df_plot,means)
      
ax.grid(which='minor')
#ax.legend(["lu", "cg","gmres", "gamg", "cg_boomeramg"])
#ax.legend(["lu", "cg","gamg"])
ax.legend(["lu", "cg","gmres","gamg"])
ax.set_xlabel("Elements")
ax.set_ylabel("Duration (s)")
#ax.legend(["lu-100","lu-200","lu-500","lu-1000","lu-1480","cg-100", "cg-200", "cg-500", "cg-1000", "cg-1480", "gamg-100", "gamg-200", "gamg-500", "gamg-1000", "gamg-1480"])

#ax2.legend(["lu", "cg", "gamg"])
ax2.legend(["lu", "cg","gmres", "gamg"])
ax2.set_xlabel("Elements")
ax2.set_ylabel("Iterations")
ax2.grid(which='minor')
plt.show()

#df_plot.info()
#df_plot.describe()
#print(df_plot)

#df_plot.plot(y="durationSolve_implicitSolver",logx=True, logy=True, marker='o')
#plt.show()
