#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Compare output data to analytical solution
# Arguments: [<show_plots (0 or 1)>] [<filenames>]
#

import sys
import numpy as np
import csv
import collections
import copy
import os
from argparse import ArgumentParser
import time
#sys.path.append('../../scripts/')
import py_reader
import json

files = ""
parser=ArgumentParser()
parser.add_argument('--path')
args=parser.parse_args()

if args.path is not None:
  path = args.path
else:
  path = './'

# get all input data in current directory
ls = os.listdir(path)

# sort files by number in file name
files = sorted(ls)

# extract the files that are npy files
solution_condition = lambda filename: ".py" in filename
solution_files = list(np.extract(np.array(list(map(solution_condition, files))), files))

print("{} files".format(len(solution_files)))

solution_files = sorted(solution_files)
solution_files_with_path = []
data = []

for solution in solution_files:
  solution_with_path = path + solution
  solution_files_with_path.append(solution_with_path)

print(solution_files_with_path)
# load py files of current group
for solution_with_path in solution_files_with_path:
  with open(solution_with_path,'rt') as f:
       dict_from_file = json.load(f)
  data.append(dict_from_file)

if len(data) == 0:
    print("no data found.")
    sys.exit(0)

####################
# 1D

with open(path + 'snapshots.csv', "wt") as csvfile:
    csvwriter = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    for dataset in data:
      ydata = py_reader.get_values(dataset, "solution", "0")
      csvwriter.writerow(ydata)
    csvfile.close()
    
sys.exit(0)
