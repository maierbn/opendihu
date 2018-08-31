#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# read csv file
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

with open('data.csv', "r") as csvfile:
    csvreader = csv.reader(csvfile, delimiter=',')
    print(type(csvreader))
    # csv_data = csvreader.astype(np.float)
    csvfile.close()
sys.exit(0)
