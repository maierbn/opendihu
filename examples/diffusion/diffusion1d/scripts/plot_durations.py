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

input_filename = "logs/log.csv"

show_plots = True
if len(sys.argv) > 1:
  show_plots = False;

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
print("df info:")
df.info()
print(df.head())


print("n rows:")
print(len(df.index))

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

plt.figure()
for solver in ["lu", "cg","gmres"]:
  df = df_all[df_all["implicitSolver_solverType"] == solver]

  means = df.groupby(['~nElements']).mean()
  print(means)

  errors = df.groupby(['~nElements']).std()
  print(errors)

  ax = means.plot(y="durationSolve_implicitSolver",logx=True, logy=True, marker='o')
  
for solver in ["gamg"]:
  #df = df_all[df_solver[df_all["implicitSolver_solverType"] == "preonly" and df_all["implicitSolver_preconditionerType"== solver]]]  
  df = df_all[df_all["implicitSolver_preconditionerType"]== solver]

  means = df.groupby(['~nElements']).mean()
  print(means)

  errors = df.groupby(['~nElements']).std()
  print(errors)

  ax = means.plot(y="durationSolve_implicitSolver",logx=True, logy=True, marker='*')
  
plt.grid(which='major')
plt.show()

