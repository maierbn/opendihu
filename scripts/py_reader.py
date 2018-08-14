#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Functions to parse opendihu *.py output files
#

import pickle, json

def get_values(data, field_variable_name, component_name):
  """
    extract the values of a single component of a field variable
    :param data: a single dict containing the data
    :param field_variable_name: the name of the field variable to consider
    :param component_name: the name of the component of the field_variable. If "0" is given, it takes the first component of the field variable.
  """
  
  for field_variable in data['data']:
    if field_variable['name'] == field_variable_name:
      if component_name == "0":
        component_name = field_variable['components'][0]['name']
      for components in field_variable['components']:
        if components['name'] == component_name:
          values = components['values']
          return values
  
def get_min_max(data, field_variable_name, component_name):
  """
    get the minimum and maximum of a field_variable and component 
    :param data: list of dicts, from multiple input files
    :param field_variable_name: the name of the field variable to consider
    :param component_name: the name of the component of the field_variable. If "0" is given, it takes the first component of the field variable.
  """
  
  min_value = None
  max_value = None
  
  # find minimum and maximum solution values
  for item in data:
    for field_variable in item['data']:
      if field_variable['name'] == field_variable_name:
        if component_name == "0":
          component_name = field_variable['components'][0]['name']
        for components in field_variable['components']:
          if components['name'] == component_name:
            values = components['values']
            item_min = min(values)
            item_max = max(values)
            
            if min_value == None or item_min < min_value:
              min_value = item_min
            if max_value == None or item_max > max_value:
              max_value = item_max
            break
        break
  return min_value, max_value
  
def load_data(filenames):
  """
    load raw data structures from opendihu *.py output files
    :param filenames: a list with the filenames
    :return: a list of dicts from the files
  """
  data = []

  # load py files
  for filename in filenames:

    # try to load file content using json, this works if it is an ascii file
    try:
      with open(filename,'r') as f:
        dict_from_file = json.load(f)
    except Exception as e:
      # try to load file contents using pickle, this works for binary files
      try: 
        with open(filename,'rb') as f:
          dict_from_file = pickle.load(f)
      except:
        
        # loading did not work either way
        dict_from_file = {}
        print("Could not parse file \"{}\"".format(filename))
        continue
    
    #print "file: {}, dict: {}".format(filename,dict_from_file)

    data.append(dict_from_file)
  return data


def get_field_variable_names(data):
  """
    get the names of the field variables contained in data
    :param data: a single dict containing the data
  """
  
  result = []
  for field_variable in data[u'data']:
    result.append(field_variable['name'])
  return result
  

def get_component_names(data, field_variable_name):
  """
    get the component names of the components of the field variable
    :param data: a single dict containing the data
    :param field_variable_name: the name of the field variable
  """
  
  result = []
  for field_variable in data['data']:
    if field_variable['name'] == field_variable_name:
      for components in field_variable['components']:
        result.append(components['name'])
  return result
  
