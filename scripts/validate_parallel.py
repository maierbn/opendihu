#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Script that parses parallel (file.0.py, file.1.py, ...) and serial (file.py) files and checks if the data matches.
#

import sys, os
import numpy as np
import py_reader    # reader utility for opendihu *.py files

all_tests_successful = True

# check if the files given in the list parallel_filenames  and the single serial_filename have the same values
def check_files(base_filename, serial_filename, parallel_filenames):
  global all_tests_successful
  #print("{}: check_files {} {}".format(base_filename, serial_filename, parallel_filenames))
  
  if serial_filename is None:
    print("{}: Serial file not found!".format(base_filename))
    return
  
  if len(parallel_filenames) == 0:
    print("{}: Parallel files not found!".format(base_filename))
    return
  
  data_parallel = py_reader.load_data(parallel_filenames)[0]
  data_serial = py_reader.load_data([serial_filename])[0]
  
  # check that non-data fields match
  field_names = ['dimension', 'nElements', 'timeStepNo', 'currentTime', 'basisFunction', 'meshType', 'basisOrder', 'onlyNodalValues']
  
  files_are_equal = True
  for field_name in field_names:
    if data_parallel[field_name] != data_serial[field_name]:
      print("{}: mismatch in '{}': serial={} != parallel={}".format(base_filename, field_name, data_serial[field_name], data_parallel[field_name]))
      files_are_equal = False
  
  # check that data is the same
  for field_variable_serial, field_variable_parallel in zip(data_serial['data'], data_parallel['data']):
    if field_variable_serial['name'] != field_variable_parallel['name']:
      print("{}: mismatch in field variable name: serial={} != parallel={}".\
        format(base_filename, field_variable_serial['name'], field_variable_parallel['name']))
      files_are_equal = False
  
    for component_serial, component_parallel in zip(field_variable_serial['components'], field_variable_parallel['components']):
      if component_serial['name'] != component_parallel['name']:
        print("{}: mismatch in component name for field variable {}: serial={} != parallel={}".\
          format(base_filename, field_variable_serial['name'], component_serial['name'], component_parallel['name']))
        files_are_equal = False
      
      for i, (value_serial, value_parallel) in enumerate(zip(component_serial['values'],component_parallel['values'])):
        diff = abs(value_serial - value_parallel)
        tolerance = 1e-2
        if diff > tolerance:
          print("{}: mismatch in \"{}\".{}[{}]: serial={} != parallel={} (|diff|={})".\
            format(base_filename, field_variable_serial['name'], component_serial['name'], i, value_serial, value_parallel, diff))
          files_are_equal = False
     
  if files_are_equal:
    print("{} ({} parallel files): match".format(base_filename, len(parallel_filenames)))
  else:
    all_tests_successful = False
  #[{u'timeStepNo': -1, u'basisFunction': u'Lagrange', 'nElements': [3, 2], u'nRanks': 2, u'meshType': u'StructuredDeformable', 
  #u'basisOrder': 1, u'currentTime': 0.0, u'onlyNodalValues': True, u'data': [{u'name': u'geometry', u'components': [{u'values': array([ 0.,  2.,  4.,  6.,  0.,  2.,  4.,  6.,  0.,  2.,  4.,  6.]), u'name': u'x'}, {u'values': array([ 0.,  0.,  0.,  0.,  2.,  2.,  2.,  2.,  4.,  4.,  4.,  4.]), u'name': u'y'}, {u'values': array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.]), u'name': u'z'}]}, {u'name': u'solution', u'components': [{u'values': array([  0.        ,   1.        ,   2.        ,   3.        ,
  #       4.24285714,   5.97142857,  10.52857143,  12.25714286,
  #       0.        ,  10.        ,  20.        ,  30.        ]), u'name': u'0'}]}, {u'name': u'rhs', u'components': [{u'values': array([  0.        ,   1.        ,   2.        ,   3.        ,
  #      -3.66666667, -11.        , -22.        , -12.83333333,
  #       0.        ,  10.        ,  20.        ,  30.        ]), u'name': u'0'}]}], u'dimension': 2}]

files = ""

if len(sys.argv) > 1:
  files = sys.argv[1:]
else:
  # get all input data in current directory
  ls = os.listdir(".")

  # sort files by number in file name
  files = sorted(ls)

# extract the files that are npy files
solution_py_condition = lambda filename: ".py" in filename
solution_py_files = list(np.extract(np.array(list(map(solution_py_condition, files))), files))

# sort files by number in file name
solution_py_files = sorted(solution_py_files)

# output file names
if len(solution_py_files) == 1:
  print("1 file: {}".format(solution_py_files[0]))
elif len(solution_py_files) <= 4:
  print("{} files: {}".format(len(solution_py_files),solution_py_files))
else:
  print("{} files: {}, {}, {}, ..., {}".format(len(solution_py_files), solution_py_files[0],solution_py_files[1],solution_py_files[2],solution_py_files[-1]))

current_timestep_parallel_filenames = []   # list of parallel filenames of the same timestep
current_timestep_serial_filename = None
current_timestep_base_filename = None        # the base filename ("file" in "file.0.py") 

for filename in solution_py_files:
  rank_no = 0
  pos_last_point = filename.rfind(".py")
  pos_2ndlast_point = filename.rfind(".", 0, pos_last_point)
  
  # if this is a parallel filename
  if pos_last_point != -1 and pos_2ndlast_point != -1:      
    base_filename = filename[0:pos_2ndlast_point]
    
    # if the currently parsed file matches the previous parsed file
    if current_timestep_base_filename != base_filename:
      if current_timestep_base_filename is not None:
        check_files(current_timestep_base_filename, current_timestep_serial_filename, current_timestep_parallel_filenames)
      
      # clear current_timestep_parallel_filenames
      current_timestep_parallel_filenames = []
      current_timestep_base_filename = base_filename
      
    current_timestep_parallel_filenames.append(filename)
  
  # if this is a serial filename
  else:
    base_filename = filename[0:pos_last_point]
    
    if current_timestep_base_filename != base_filename:
      if current_timestep_base_filename is not None:
        check_files(current_timestep_base_filename, current_timestep_serial_filename, current_timestep_parallel_filenames)
      
      # clear current_timestep_parallel_filenames
      current_timestep_parallel_filenames = []
      current_timestep_base_filename = base_filename
      
    current_timestep_serial_filename = filename

check_files(current_timestep_base_filename, current_timestep_serial_filename, current_timestep_parallel_filenames)

if all_tests_successful:
  sys.exit(0)
else:
  sys.exit("Files are not equal.")
