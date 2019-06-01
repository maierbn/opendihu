#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Function that compares two data sets of files and computes the relative error between the solutions.
#

import sys, os
import numpy as np
import py_reader    # reader utility for opendihu *.py files

def compute_error_filenames(reference_filenames, filenames_to_test, field_variable_name="solution"):
  
  if reference_filenames is None or len(reference_filenames) == 0:
    print("No reference files given")
    return None
    
  if filenames_to_test is None or len(filenames_to_test) == 0: 
    print("No files to test given")
    return None
    
  data_sets_reference = py_reader.load_data(reference_filenames)
  data_sets_to_test = py_reader.load_data(filenames_to_test)
    
  if len(data_sets_reference) != len(data_sets_to_test):
    print("Number of timesteps does not match! Reference data set: {}, Test data set: {}".format(data_sets_reference, data_sets_to_test))
    
  data_reference = data_sets_reference[-1]
  data_test = data_sets_to_test[-1]
  
  # parse field variables and extract values
  values_reference = None
  values_test = None
  
  data_reference
  
  # parse reference field variable
  for field_variable in data_reference["data"]:
    if field_variable["name"] == field_variable_name:
      values_reference = field_variable["components"][0]["values"]
  
  # parse test field variable
  for field_variable in data_test["data"]:
    if field_variable["name"] == field_variable_name:
      values_test = field_variable["components"][0]["values"]
  
  if values_reference is None or values_test is None:
    return None
  
  # compute error
  absolute_error = np.linalg.norm(values_reference - values_test)
  print("values_reference: ",values_reference)
  print("values_test: ",values_test)
  print("difference: ",values_reference - values_test)
  
  try:
    relative_error = np.linalg.norm((values_reference - values_test) / values_reference)
  except:
    relative_error = np.linalg.norm((values_reference - values_test) / values_test)
  
  print("Absoluter Fehler: ", absolute_error)

  relative_error /= len(values_reference)

  print("Relativer Fehler: ",relative_error)
  
  return relative_error
  
def get_matching_py_filenames(part_of_filename):
  """
  given a file base like "a/b/c/file_" return all matching *.py files as a list, e.g. ["a/b/c/file_0.py", "a/b/c/file_1.py"]
  """  
  # get directory
  path = os.path.dirname(part_of_filename)
  if path == "":
    path = "."
  part_of_filename = os.path.basename(part_of_filename)
    
  # get all files in path, sorted
  files = sorted(os.listdir(path))
  
  condition = lambda filename: part_of_filename in filename and ".py" in filename
  matching_files = list(np.extract(np.array(list(map(condition, files))), files))
  
  result = []
  for matching_file in matching_files:
    result.append(os.path.join(path,matching_file))
    
  return result

def compute_error(reference_filenames_base, files_to_test_base, field_variable_name="solution"):
  """
  return the relative error, divided by the number of entries between two opendihu python output files (or sets of files, like file_0.0.py, file_0.1.py)
  """  
  reference_filenames = get_matching_py_filenames(reference_filenames_base)
  test_filenames = get_matching_py_filenames(files_to_test_base)
  print("reference_filenames_base: {}, reference_filenames: {}".format(reference_filenames_base,reference_filenames))
  print("files_to_test_base: {}, test_filenames: {}".format(files_to_test_base,test_filenames))
  
  return compute_error_filenames(reference_filenames, test_filenames, field_variable_name)

if __name__ == '__main__':
  compute_error(sys.argv[1], sys.argv[2], "solution")

