#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# pandas utility to load data frames from csv files like they are produced by opendihu. 
# use df = pandas_utility.load_df(input_filename) to get the data from a csv file
# then directly call pandas_utility.print_table() and pandas_utility.plot_weak_scaling() if needed
# 
# If this script is executed directly, it plots the given filename as weak scaling data set.
# usage: ./pandas_utility.py <input_filename log.csv>
# 


# good read about how pandas is used: https://towardsdatascience.com/pandas-equivalent-of-10-useful-sql-queries-f79428e60bd9

import sys, os
import subprocess
import datetime
import time
import numpy as np
import pandas as pd
import json
from cycler import cycler
# load matplotlib
import matplotlib
import matplotlib.pyplot as plt

def merge_dicts(x, y):
  z = x.copy()   # start with x's keys and values
  z.update(y)    # modifies z with y's keys and values & returns None
  return z

def remove_duplicates(seq):
  seen = set()
  seen_add = seen.add
  return [x for x in seq if not (x in seen or seen_add(x))]

def determine_column_names(line):
  column_names = list(line.split(";"))

  # rename "n" columns 
  for i, column_name in enumerate(column_names):
    if column_name == "n":
      column_names[i] = "{}_n".format(column_names[i-1])

  while "" in column_names:
    column_names.remove("")
  while "\n" in column_names:
    column_names.remove('\n')
  while "0" in column_names:
    column_names.remove('0')

  # remove duplicates
  seen = set()
  column_names2 = []
  for x in list(column_names):
    if x not in seen:
      seen.add(x)
      column_names2.append(x)
    else:
      print("  Note: column \"{}\" appears multiple times".format(x))
      while x in seen:
        x = "{}_2".format(x)
      seen.add(x)
      column_names2.append(x)

  column_names = column_names2
  return column_names
  n_columns = len(column_names)

def load_df(input_filename):
  """
  Create a data frame from the csv file given as input. This is a robust version of the pandas pd.read_cvs() function that handles missing/invalid lines etc.
  It allows different data sets to have different columns, if there is a new column header each time (as there is in opendihu generated csv files).
  The timestamp is converted to datetime format, also a new column duration_init is computed.
  """
  
  # if the pickle file does not yet exist, try to remove the .pickle suffix
  if not os.path.isfile(input_filename) and ".pickle" in input_filename:
    input_filename = input_filename[0:input_filename.index(".pickle")]
  
  if ".pickle" in input_filename:
    
    print("Read pickle file \"{}\"".format(input_filename))
    df = pd.read_pickle(input_filename)
    
  elif ".json" in input_filename:
    
    
    with open(input_filename) as f:
      for i, l in enumerate(f):
        pass
    n_lines = i+1
    
    print("Parsing json in file \"{}\", {} lines".format(input_filename, n_lines))
    
    df_lines = []
    with open(input_filename) as f:
      for i,line in enumerate(f):
        if "{" not in line:
          continue
        
        if line[0] == "#":
          continue
          
        # parse json code in line
        python_dict = json.loads(line)
        new_python_dict = {}
        
        for (key,value) in python_dict.items():
          try:
            new_python_dict[key] = [float(value)]
          except:
            new_python_dict[key] = [value]
        #python_dict = {key: [value] for (key,value) in python_dict.items()}

        # create data frame of only this line
        df = pd.DataFrame.from_dict(new_python_dict)
        
        df_lines.append(df)
        
        # print progress
        if i % 100 == 0:
          sys.stdout.write("{:.2f}%\b\b\b\b\b\b".format(100*i/n_lines))
        
    print("File \"{}\" contains {} lines, creating data frame.".format(input_filename, len(df_lines)))
        
    # concat all lines
    df = pd.concat(df_lines)
  else:
    
    # determine columns to load from the log file
    with open(input_filename) as f:
      
      # separate lines of file in blocks for a single scenario
      scenario_blocks = []
      current_lines = []
      column_names = None
      for line in f:
        if "~" in line:
          if current_lines != []:
            scenario_blocks.append(current_lines)
          current_lines = []
        current_lines.append(line)
      scenario_blocks.append(current_lines)
        
      print("File \"{}\" contains data from {} runs".format(input_filename,len(scenario_blocks)))
        
      # loop over scenario blocks
      df_blocks = []
      for scenario_block in scenario_blocks:
        if "~" in scenario_block[0] or True:
          column_names = determine_column_names(scenario_block[0])
          n_columns = len(column_names)
          
          #print("parse block with {} columns".format(n_columns))
          with open(".pd","w") as fout:
            for line in scenario_block:
              fout.write(line)
                                  
        # load data frame
        try:
          df_block = pd.read_csv(".pd", sep=';', error_bad_lines=False, warn_bad_lines=True, comment="#", header=None, names=column_names, usecols=range(n_columns), mangle_dupe_cols=True)
          
          # disable console output about loading blocks
          if False:
            if df_block.shape[0] == 0:
              print("skip empty block")
            elif "scenarioName" in df_block:
              print("load block of {} rows for scenario {}".format(df_block.shape[0], df_block.iloc[0]["scenarioName"]))
            else:
              print("load block of {} rows".format(df_block.shape[0]))
          
          df_blocks.append(df_block)
        except Exception as e:
          print("could not load block: {}".format(str(e)))
      
    # concat the blocks of all scenarios
    df = pd.concat(df_blocks)
    
  # preprocessing
  # --------------
    
  # parse timestamp
  if "timestamp" in df:
    df["# timestamp"] = df["timestamp"]
  
  if "# timestamp" in df:
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

  if "duration_total_advanceTimeSpan1" in df and "duration_0D_fast" in df and "duration_1D_fast" in df:
    df['duration_transfer_01D'] = df['duration_total_advanceTimeSpan1'] - df['duration_0D_fast'] - df['duration_1D_fast']
  elif "duration_total_advanceTimeSpan1" in df and "duration_0D" in df and "duration_1D" in df:
    df['duration_transfer_01D'] = df['duration_total_advanceTimeSpan1'] - df['duration_0D'] - df['duration_1D']
 
  # save df as pickle if it was no pickle
  if ".pickle" not in input_filename:
    output_filename = "{}.pickle".format(input_filename)
    df.to_pickle(output_filename)
    print("Saved to \"{}\".".format(output_filename))
 
  n_columns = len(df.columns)
  print("Loaded {} rows and {} colums.".format(len(df.index), n_columns))
  
  return df

def compute_means_stddev(df, column_names,minmax=False):
  """
  compute items for computation of means and standard deviation
  usage is as follows:
  
  mean_items,stddev_items = compute_means_stddev(df, columns_to_plot)
    
  # aggregate data
  means = df.groupby(group_by_columns).agg(mean_items)
  errors = df.groupby(group_by_columns).agg(stddev_items)
  """
  
  column_names = list(set(column_names))
  
  # create auxiliary columns that will be computed
  if not "memoryResidentSet" in df:
    df["memoryResidentSet"] = 0
  df["n"] = 0

  # remove column that are not present in the df
  for column in list(column_names):
    if column not in df:
      print("Note: remove invalid column {} from column_names".format(column))
      column_names.remove(column)
  
  # define items to be printed, the columns "n" and "memoryResidentSet" need to be already present in the df
  means = merge_dicts(
    {column_name: lambda v: (np.mean(v) if v.dtype == np.float64 else str(v.iloc[0]) ) for column_name in column_names},
    {'n': np.size, "memoryResidentSet": lambda v: "{:.3f} GB".format(np.mean(v)/(1024.**3))}
  )
  
  stddevs = merge_dicts(
    {column_name: lambda v: (np.std(v) if v.dtype == np.float64 else str(v.iloc[0]) ) for column_name in column_names},
    {'n': np.size, "memoryResidentSet": lambda v: "{:.3f} GB".format(np.std(v)/(1024.**3))}
  )

  if minmax:
    mins = merge_dicts(
      {column_name: lambda v: (np.min(v) if v.dtype == np.float64 else str(v.iloc[0]) ) for column_name in column_names},
      {'n': np.size, "memoryResidentSet": lambda v: "{:.3f} GB".format(np.min(v)/(1024.**3))}
    )

    maxs = merge_dicts(
      {column_name: lambda v: (np.max(v) if v.dtype == np.float64 else str(v.iloc[0]) ) for column_name in column_names},
      {'n': np.size, "memoryResidentSet": lambda v: "{:.3f} GB".format(np.max(v)/(1024.**3))}
    )
    
    return means, stddevs, mins, maxs
  return means, stddevs
  
def print_table(df, title, columns_to_print_longnames, column_shortnames, group_by_columns=['scenarioName','nRanks']):
  """
  Print values of a df to console. Compute mean values.
  
  :param df:    the pandas dataframe to use
  :param title: title that is printed before the table
  :param columns_to_print_longnames: a list of names of columns that should be included in the table
  :param column_shortnames: a dict that maps long names to short names. The short names will be the caption of the columns, if given. 
                   It makes sense to specify short names for long column names to get more information in the table for the given maximum width.
  :param group_by_columns: the columns that define one data point
  """
  
  # compute mean values
  mean_items,_ = compute_means_stddev(df, columns_to_print_longnames)

  # aggregate data
  agg_df = df.groupby(group_by_columns).agg(mean_items).rename(columns = column_shortnames)

  available_columns_shortnames = list(agg_df.columns)
  columns_to_print_shortnames = [column_shortnames[long_name] if long_name in column_shortnames else long_name for long_name in columns_to_print_longnames]

  #print("available colums: {}".format(available_columns_shortnames))
  #print("columns_to_print_longnames: {}".format(columns_to_print_longnames))
  #print("columns_to_print_shortnames: {}".format(columns_to_print_shortnames))
  
  # remove columns that are not present in columns_to_print_longnames or available_columns_shortnames
  for shortname,longname in list(zip(columns_to_print_shortnames,columns_to_print_longnames)):
    if longname not in columns_to_print_longnames or shortname not in available_columns_shortnames:
      columns_to_print_shortnames.remove(shortname)

  # print values as a table
  print("-"*120)
  print(title)
  print(agg_df[columns_to_print_shortnames])
  print("-"*120)

def plot_weak_scaling(df, title, columns_to_plot, group_by_columns=['nRanks'], plot_labels=None):
  """
  Make a weak scaling plot out of the df, write it to a png and pdf file in `out` directory.
  You have to call plt.show() to show the plotting window afterwards.
  Mean and standard deviation values are computed.
  
  :param df:    the pandas dataframe to use
  :param title: title and output file name of the plot
  :param columns_to_plt: a list of names of columns that should be one line each
  :param group_by_columns: the columns that define one data point
  :param plot_labels: a list of strings, custom labels to be shown in the legend
  """
  
  # compute mean and stddev values
  mean_items,stddev_items,min_items,max_items = compute_means_stddev(df, columns_to_plot, minmax=True)
  if len(mean_items) == 0:
    return
    
  # aggregate data
  means = df.groupby(group_by_columns).agg(mean_items)
  errors = df.groupby(group_by_columns).agg(stddev_items)
  mins = df.groupby(group_by_columns).agg(min_items)
  maxs = df.groupby(group_by_columns).agg(max_items)
  
  means = means.interpolate(method='linear', limit_area='inside')

  # remove not existent column names from columns_to_plot
  available_columns_longnames = list(means.columns)
  
  for column_name in list(columns_to_plot):
    if column_name not in available_columns_longnames:
      print("Note: remove invalid column {} from columns_to_plot".format(column_name))
      columns_to_plot.remove(column_name)
    
  # define colors and linestyles, do this befor this function
  #plt.rc('axes', prop_cycle=(cycler('color', ['y', 'r', (0.1,0.1,0.6), (0.5,0.5,0.5), (1.0,0.7,0.2)]) +
  #                           cycler('linestyle', [ '-', '-', '-', '--','-'])))
    
  # plot values
  print("<<<<<<<means")
  print(means)
  print("<<<<<<<stddev")
  print(errors)
  print("<<<<<<<min")
  print(mins)
  print("<<<<<<<max")
  print(maxs)
  
  # create single long table for pgf plots
  df_pgf = pd.concat([means.drop(columns=['memoryResidentSet']).add_prefix('mean_'), mins.drop(columns=['memoryResidentSet']).add_prefix('min_'), maxs.drop(columns=['memoryResidentSet']).add_prefix('max_')], axis=1)
  with open("out.pgf", "w") as f:
    f.write("# {}\n".format(title))
    f.write("# created by pandas_utility.py at {}\n".format(datetime.datetime.now().strftime('%d.%m.%Y %H:%M:%S')))
    f.write(str(df_pgf))
  print("Write file \"out.pgf\".")
   
  # plot with std deviation
  #ax = means.plot(figsize=(13,7), y=columns_to_plot, title = title, logx=True, logy=True, yerr=errors, marker='o'))
  
  # plot with min max
  lower_error_bars = pd.DataFrame()
  upper_error_bars = pd.DataFrame()
  print(columns_to_plot)
  for column in columns_to_plot:
    lower_error_bars[column] = means[column]-mins[column]
    upper_error_bars[column] = maxs[column]-means[column]
  
  l = lower_error_bars.to_numpy()[:, np.newaxis, :]
  h = upper_error_bars.to_numpy()[:, np.newaxis, :]
  
  yerr = np.concatenate((l,h),axis=1)
  yerr = np.moveaxis(yerr, [0,2], [2,0])
  
  ax = means.plot(figsize=(13,7), y=columns_to_plot, title = title, logx=True, logy=True, yerr=yerr, marker='o')
  
  # set axis labels
  rank_nos = sorted(list(set(df["nRanks"])))
  print("rankNos on x axis: {}".format(rank_nos))
  
  rank_no_labels = []
  for i in rank_nos:
    
    if i == 7776:
      rank_no_labels.append("\n7776")
    elif i == 27744:
      rank_no_labels.append("\n27744")
    else:
      rank_no_labels.append(str((int)(i)))
  
  ax = plt.gca()
  ax.set_xticks(rank_nos)
  ax.set_xticklabels(rank_no_labels)
  ax.grid(which='major')
  ax.set_xlabel('Number of processes')
  ax.set_ylabel('Runtime [s]')
  ax.set_title(title)
  if plot_labels is not None:
    plt.subplots_adjust(right=0.66, top=0.94, bottom=0.18)
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., frameon=False, labels=plot_labels)  
  else:
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., frameon=False)
  
  plt.tight_layout()
  filename1 = "out/{}.pdf".format(title.replace(" ", ""))
  filename2 = "out/{}.png".format(title.replace(" ", ""))
  subprocess.check_output(["mkdir", "-p", "out"])
  
  print("write files {} and {}".format(filename1,filename2))
  plt.savefig(filename1)
  plt.savefig(filename2)
 

if __name__ == "__main__":
  # execute only if run as a script
  
  input_filename = "logs/log.csv"
  if len(sys.argv) == 1:
    print("usage: ./pandas_utility.py <input_filename log.csv> [-l] [-p] [-h] [<field names to display>]\n  -l: list field names\n  -p: create weak scaling plot\n  -h: show help\n")
    quit()
  else:
    input_filename = sys.argv[1]
    
  list_columns = False
  create_plot = False
  custom_columns_to_print = []
  
  if len(sys.argv) >= 3:
    for arg in sys.argv[2:]:
      if arg == "-h" or arg == "--help":
        print("usage: ./pandas_utility.py <input_filename log.csv> [-l] [-p] [-h] [<field names to display>]\n  -l: list field names\n  -p: create weak scaling plot\n  -h: show help\n")
        quit()
      if arg == "-l":
        list_columns = True
        continue
      elif arg == "-p":
        create_plot = True
        continue
      custom_columns_to_print.append(arg)
  
  if not custom_columns_to_print:
    print("Note: You can specify which columns to display by:\n  ./pandas_utility.py <input_filename log.csv> [-l] [-p] [-h]  [<field names to display>]")
    print('Use the "-l" option to learn all possible columns.\n')
  
  
  # ------------------------------------------------
  # define shortnames for the table, each line is
  #  long_name : short_name
  column_shortnames = {
    "totalUsertime": "user",
    "duration_total": "total comp.",
    "duration_0D": "0D",
    "duration_1D": "1D",
    #"duration_init": "duration_init",
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
    "nIterations_multidomainLinearSolver": "niter",
  }

  # define columns for table and plot (long names)
  columns_to_print = [
    "meta_partitioning", 
    "totalUsertime", 
    "duration_total",
    "duration_init",  
    "durationOnlyWrite",
    "duration_0D", 
    "duration_1D", 
    "duration_bidomain", 
    "duration_multidomain", 
    "duration_mechanics", 
    "nIterations_multidomainLinearSolver",
    "memoryResidentSet", 
    "n"
  ]
  
  if custom_columns_to_print:
    columns_to_print = custom_columns_to_print
  
  columns_to_plot = [
    "totalUsertime", "duration_total", "duration_bidomain", "duration_multidomain", "duration_mechanics", 
    "duration_0D", "duration_1D"
  ]

  # load data frame
  df = load_df(input_filename)
  
  if list_columns:
    print("all columns:")
    for column in df.columns:
      print("  "+column)
  
  title = ""
  columns_to_plot = ["totalUsertime", "duration_total", "duration_bidomain", "duration_0D", "duration_1D", "duration_transfer_01D"]
  plot_labels = ["user time", "total computation", "3D model", "0D model", "1D model", "communication 0D,1D"]

  # print table
  print_table(df, title, columns_to_print, column_shortnames)

  # plot weak scaling
  if create_plot:
    title = "Weak scaling"
    plot_weak_scaling(df, title, columns_to_plot, ['nRanks'], plot_labels)
  
    plt.show()
