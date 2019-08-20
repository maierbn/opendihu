import sys

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


# example
# serie_lu = pd.Series({"100": 0.0126655, "200": 0.0115678, "500":0.0116677, "1000": 0.0128555, "1480": 0.0124516}, name = "lU")

if(len(sys.argv)) >= 2:
  input_filename = sys.argv[1]
else:
  input_filename = "logs/log.csv"
print("filename: {}".format(input_filename))

def load_df(input_filename):

  # determine number of columns
  with open(input_filename) as f:
    line = f.readline()
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

df_duration_mean = pd.DataFrame();

#for solver in ["lu", "cg","gmres", "gamg", "cg_boomeramg"]:
for solver in ["lu"]:
  df_name = "df_" + solver
  print("df: {}".format(df_name))
  
  df_solver = df[df["scenarioName"] == ("serial_diffusion_" + solver)]
  #df_solver.info()
  
  #df_solver_duration = df_solver[df_solver["durationSolve_implicitSolver"]]
  #df_solver_duration.info()
  #df_solver_duration.describe()
  
  df_solver_means = df_solver.groupby(['~nElements']).mean()
  print(df_solver_means)
  
  
  #duration_error = df_solver["durationSolve_implicitSolver"].groupby(['~nElements']).std() 
  
  #df_duration_mean.append(duration_mean)

df_solver_means.plot(y = "totalUsertime", logx=True, logy=True, marker='o')

plt.grid(which='major')
plt.show()
  
