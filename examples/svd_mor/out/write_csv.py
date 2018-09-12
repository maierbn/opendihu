#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Write csv file of solution for all timesteps
# Arguments: [<show_plots (0 or 1)>] [<filenames>]
#

import sys
import numpy as np
import collections
import copy
from sets import Set
import os
import time
import datetime
import scipy.integrate

import py_reader
import csv


files = ""
# get all input data in current directory
ls = os.listdir(".")

# sort files by number in file name
files = sorted(ls)

# extract the files that are npy files
solution_condition = lambda filename: ".py" in filename
solution_files = list(np.extract(map(solution_condition, files), files))

print("{} files".format(len(solution_files)))

data = py_reader.load_data(solution_files)

if len(data) == 0:
    print("no data found.")
    sys.exit(0)

dimension = data[0]['dimension']


with open('data.csv', "a") as csvfile:
    csvwriter = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

    for dataset in data:
	timestep_no = dataset['timeStepNo']
	t = dataset['currentTime']
	ydata = py_reader.get_values(dataset, "solution", "0")
	csvwriter.writerow(ydata)
    csvfile.close()

sys.exit(0)
