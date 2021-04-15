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
import pandas_utility

input_filename = "logs/log.csv"
list_columns = False
show_plots = False

# parse arguments
for arg in sys.argv[1:]:
  if ".csv" in arg:
    input_filename = arg
  else:
    if "l" in arg:
      list_columns = True
    if "p" in arg:
      show_plots = True

if len(sys.argv) == 1:
  print("usage: {} [-l] [-p] [<input_filename>]".format(sys.argv[0]))
  print("   -l: list column names")
  print("   -p: show plot window (else create pdf plot)")
  print("   <input_filename> log file, default: logs/log.csv")
  print("")


# load matplotlib
import matplotlib
if not matplotlib.is_interactive() or not show_plots:
  matplotlib.use('Agg')

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# define global plotting parameters
plt.rcParams.update({'font.size': 16})
plt.rcParams['lines.linewidth'] = 2

def remove_duplicates(seq):
  seen = set()
  seen_add = seen.add
  return [x for x in seq if not (x in seen or seen_add(x))]

# determine columns to load from the log file
with open(input_filename) as f:
  line = f.readline()
  if "~nDofs" in line:
    pos = line.find("~nDofs")
    line = line[0:pos]
  column_names = list(line.split(";"))

  # rename "n" columns 
  for i, column_name in enumerate(column_names):
    if column_name == "n":
      column_names[i] = "{}_n".format(column_names[i-1])

  while "" in column_names:
    column_names.remove("")

  seen = set()
  column_names2 = []
  for x in list(column_names):
    if x not in seen:
      seen.add(x)
      column_names2.append(x)
    else:
      print("Note: column \"{}\" appears multiple times".format(x))
      while x in seen:
        x = "{}_2".format(x)
      seen.add(x)
      column_names2.append(x)

  column_names = column_names2
  n_columns = len(column_names)
  
if list_columns:
	print("File {} contains {} colums: {}".format(input_filename, n_columns, column_names))

# load data frame
#df = pd.read_csv(input_filename, sep=';', error_bad_lines=False, warn_bad_lines=True, comment="#", header=None, names=column_names, usecols=range(n_columns), mangle_dupe_cols=True)
df = pandas_utility.load_df(input_filename)

# filter data
#df = df.loc[df['endTime'] == 1]

if not list_columns:
  print("File {} contains {} rows and {} colums.".format(input_filename, len(df.index), n_columns))

# parse timestamp
df['# timestamp'] =  pd.to_datetime(df['# timestamp'])

# compute new field for initialization time
if 'durationParaview3DInit' not in df:
    df['durationParaview3DInit'] = 0
if 'durationParaview1DInit' not in df:
    df['durationParaview1DInit'] = 0
if 'durationParaview1DWrite' not in df:
    df['durationParaview1DWrite'] = 0
    
if 'durationParaview3DWrite' not in df:
  if 'durationParaviewOutput' in df:
    df['durationOnlyWrite'] = df['durationParaviewOutput']
  else:
    df['durationOnlyWrite'] = df['durationWriteOutput']
else:
  df['durationOnlyWrite'] = df['durationParaview3DWrite'] + df['durationParaview1DWrite']

try:
  df['duration_init'] = df['totalUsertime'] - df['duration_total'] + df['durationParaview3DInit'] + df['durationParaview1DInit']
except:
  df['duration_init'] = 0

# Info about the data structure
#print("df info:")
#df.info()
#print(df.head())

# set options for console display
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
pd.set_option('display.max_colwidth', 100)

def merge_dicts(x, y):
  z = x.copy()   # start with x's keys and values
  z.update(y)    # modifies z with y's keys and values & returns None
  return z

def plot(df, title, columns_to_plot, mylabels=None):
  means = df.groupby(['nRanks']).mean()
  if len(means) == 0:
    return
  
  #linestyle_cycler = cycler.cycler('linestyle',['-','--',':','-.'])
  # (cycler.cycler('color', ["k",(0.3,0.3,0.7),(0.7,0.7,1.0), "r", "y"])+cycler.cycler('linestyle', ['-', '--', ':', '-', '-'])))
  plt.rc('axes', prop_cycle=(cycler('color', ['k', 'b', 'r', 'y']) +
                             cycler('linestyle', ['-', '--', '-.', '-'])))
  #plt.rc('axes', prop_cycle=("cycler('color', 'rgb') + cycler('linestyle',  ['-', '-', ':'])"))
    
  errors = df.groupby(['nRanks']).std()
  ax = means.plot(figsize=(10,7), y=columns_to_plot, title=title, logx=True, logy=True, yerr=errors, marker='o')
  
  rank_nos = sorted(list(set(df["nRanks"])))
  ax = plt.gca()
  ax.set_xticks(rank_nos)
  ax.set_xticklabels(rank_nos)
  ax.grid(which='major')
  ax.set_xlabel('number of processes')
  ax.set_ylabel('runtime [s]')
  ax.set_title(title)
  if mylabels is not None:
    ax.legend(labels=mylabels)
  
  if show_plots:
    plt.show()
  else:
    plt.tight_layout()
    plot_filename = "{}.pdf".format(title.replace(" ", "").replace("/", ""))
    plt.savefig(plot_filename)
    print("Created \"{}\".".format(plot_filename))
 
def output(df, title, columns_to_print, columns_to_plot, plot_labels=None):
  """
  print values to console and produce plot
  """
  columns_to_extract = list(set(columns_to_plot + columns_to_print + ["totalUsertime", "durationReadGeometry", "durationSetStiffnessMatrix", 
    "durationOnlyWrite", "durationAssembleBoundaryConditions", "durationInitCellml", "durationComputeMappingBetweenMeshes",
     "durationMap"]))

  # create auxiliary columns that will be computed
  if not "memoryResidentSet" in df:
    df["memoryResidentSet"] = 0
  df["n"] = 0

  # remove column that are not present in the df
  for column in list(columns_to_extract):
    if column not in df:
      print("Note: remove invalid column {} from columns_to_extract".format(column))
      columns_to_extract.remove(column)
  
  for column in list(columns_to_print):
    if column not in df:
      print("Note: remove invalid column {} from columns_to_print".format(column))
      columns_to_print.remove(column)
  
  for column in list(columns_to_plot):
    if column not in df:
      print("Note: remove column {} from columns_to_plot".format(column))
      columns_to_plot.remove(column)

  # select the captions for the table
  table_shortnames = [column_shortnames[long_name] if long_name in column_shortnames else long_name for long_name in columns_to_print]
  
  # define items to be printed, the columns "n" and "memoryResidentSet" need to be already present in the df
  items = merge_dicts(
    {column_name: lambda v: (np.mean(v) if v.dtype == np.float64 else str(v.iloc[0]) ) for column_name in columns_to_extract},
    {'n': np.size, "memoryResidentSet": lambda v: "{:.3f} GB".format(np.mean(v)/(1024.**3))}
  )

  print("-"*120)
  print(title)
  print(df.groupby(['scenarioName','nRanks']).agg(items).rename(columns = column_shortnames)[table_shortnames])
  print("-"*120)

  # create plot
  plot(df, title, columns_to_plot, plot_labels)

# ------------------------------------------------
# define shortnames for the table, each line is
#  long_name : short_name
column_shortnames = {
  "totalUsertime": "user",
  "duration_total": "total comp.",
  "duration_0D": "0D",
  "duration_1D": "1D",
  "duration_init": "duration_init",
  "duration_bidomain": "bidomain",
  "duration_mechanics": "mechanics",
  "duration_multidomain": "multidomain",
  "durationAssembleBoundaryConditions": "initBC",
  "durationSetStiffnessMatrix": "stiffness",
  "durationComputeMappingBetweenMeshes": "compMap",
  "durationMap": "map",
  "durationReadGeometry": "read",
  "durationOnlyWrite": "write",
  "durationInitCellml": "initCell",
  "memoryResidentSet": "mem",
  "meta_partitioning": "subdomains",
  "nonlinear_durationSolve": "nonlin",
}

# define columns for table and plot (long names)
columns_to_print = ["meta_partitioning", "totalUsertime", "duration_total", "duration_0D", "duration_1D", "duration_bidomain", "duration_multidomain", "duration_mechanics", "nonlinear_durationSolve", "duration_init", "durationOnlyWrite", "memoryResidentSet", "n"]
columns_to_plot = ["duration_total", "duration_init", "durationOnlyWrite", "duration_0D", "duration_1D", "duration_bidomain", "duration_multidomain", "duration_mechanics"]

plot_labels = ["total", "initialization", "write VTK files", "0D model", "1D model", "3D model"]

title = input_filename
output(df, title, columns_to_print, columns_to_plot, plot_labels)

