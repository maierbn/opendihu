#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Functions to parse opendihu *.py output files
#

import pickle, json
import copy
import numpy as np

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
    Load raw data structures from opendihu *.py output files.
    If there are multiple files from different processes with appropriate ending, merge the data in these files.
    :param filenames: a list with the filenames
    :return: a list of dicts from the files
  """
  loaded_data = []

  # group parallel files
  grouped_filenames = []  # each item is [<base_filename>, [<filename0>, <filename1>, ...]]
  for filename in filenames:
    if filename is None:
      continue
    rank_no = 0
    pos_last_point = filename.rfind(".py")
    pos_2ndlast_point = filename.rfind(".", 0, pos_last_point)
    
    bucket_found = False
    base_filename = filename
    if pos_last_point != -1 and pos_2ndlast_point != -1:      
      base_filename = filename[0:pos_2ndlast_point]
      rank_str = filename[pos_2ndlast_point+1:pos_last_point]
      #print("filename: {}, base_filename: {}, rank_str: {}".format(filename, base_filename, rank_str))
      #rank_no = (int)(rank_str)
    
      # find bucket in grouped_filenames with matching base_filename
      for i in range(len(grouped_filenames)):
        if grouped_filenames[i][0] == base_filename:
          grouped_filenames[i][1].append(filename)
          bucket_found = True
          break
    if not bucket_found:
      grouped_filenames.append([base_filename, [filename]])
    
  # loop over groups of filenames, one group is the files from different rank but from same scenario/time step
  for file_group in grouped_filenames:
    filenames = file_group[1]
      
    group_data = []
      
    # load py files of current group
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
          dict_from_file = None
          print("Could not parse file \"{}\"".format(filename))
          continue
      
      #print("file: {}, dict: {}".format(filename,dict_from_file))
      if dict_from_file is not None:
        group_data.append(dict_from_file)
   
    if len(group_data) == 0:
      continue
      
    # unstructured meshes are not parallel and have only a single number of elements, "nElements"
    if group_data[0]["meshType"] == "UnstructuredDeformable":
      merged_data = group_data[0]
      
    # structured meshes can be parallel and have global and local numbers of elements for each coordinate direction, "nElementsGlobal", "nElementsLocal"
    else:
     
      # merge group data
      merged_data = copy.deepcopy(group_data[0])
      del merged_data['beginNodeGlobalNatural']
      del merged_data['hasFullNumberOfNodes']
      del merged_data['nElementsLocal']
      del merged_data['ownRankNo']
      n_ranks = merged_data['nRanks']
      dimension = merged_data['dimension']

      # determine number of nodes
      average_n_nodes_per_element_1d = 1
      if merged_data['basisFunction'] == 'Lagrange':
        average_n_nodes_per_element_1d = merged_data['basisOrder']
      elif merged_data['basisFunction'] == 'Hermite':
        average_n_nodes_per_element_1d = 1

      # determine number of dofs per node
      n_dofs_per_node = 1
      if merged_data['basisFunction'] == 'Hermite' and not merged_data['onlyNodalValues']:
        n_dofs_per_node = 2**dimension

      # determine global size of arrays
      n_dofs_global = []
      n_dofs_total = 1
      for i in range(dimension):
        n_nodes_global_direction = merged_data['nElementsGlobal'][i]*average_n_nodes_per_element_1d + 1
        n_dofs_global_direction = n_nodes_global_direction*n_dofs_per_node
        n_dofs_global.append(n_dofs_global_direction)
        n_dofs_total *= n_dofs_global_direction
    
      #print("n_dofs_global: ",n_dofs_global)
    
      merged_data['nElements'] = merged_data['nElementsGlobal']
      del merged_data['nElementsGlobal']

      # loop over data fields
      for field_variable_index,field_variable in enumerate(merged_data['data']):
        for component_index,component in enumerate(field_variable['components']):
          
          # resize array
          merged_data['data'][field_variable_index]['components'][component_index]['values'] = np.zeros(n_dofs_total)

          # loop over ranks / input files
          for data in group_data:
            
            # compute local size
            n_nodes_local = []
            n_dofs_local = []
            for i in range(dimension):
              n_nodes_local_direction = data['nElementsLocal'][i]*average_n_nodes_per_element_1d
              if data['hasFullNumberOfNodes'][i]:
                n_nodes_local_direction += 1
              n_nodes_local.append(n_nodes_local_direction)
              n_dofs_local.append(n_nodes_local_direction*n_dofs_per_node)
              
            #print("")
            #print(field_variable["name"],component["name"])
            #print("data: ",data)
            #print("rank {}, n_dofs_local: {}, data: {}".format(data["ownRankNo"], n_dofs_local,data['data'][field_variable_index]['components'][component_index]['values']) )
            
            # set local portion in global array
            indices_begin = []
            indices_end = []
            for i in range(dimension):
              begin_node_global = data['beginNodeGlobalNatural'][i]
              indices_begin.append(begin_node_global*n_dofs_per_node)
              indices_end.append(begin_node_global*n_dofs_per_node+n_dofs_local[i])
            
            #print("indices {} - {}".format(indices_begin, indices_end))
            
            if dimension == 1:
              for x in range(indices_begin[0],indices_end[0]):
                index_in = x-indices_begin[0]
                index_result = x
                
                #print("x: {}, index_in: {}, index_result: {}".format(x, index_in, index_result))
                
                merged_data['data'][field_variable_index]['components'][component_index]['values'][index_result] \
                  = data['data'][field_variable_index]['components'][component_index]['values'][index_in]
            
            elif dimension == 2:
              for y in range(indices_begin[1],indices_end[1]):
                for x in range(indices_begin[0],indices_end[0]):
                  index_in = (y-indices_begin[1])*n_dofs_local[0] + (x-indices_begin[0])
                  index_result = y*n_dofs_global[0] + x
                  
                  #print("x: {}, y: {}, index_in: {}, index_result: {}".format(x, y, index_in, index_result))
                  
                  merged_data['data'][field_variable_index]['components'][component_index]['values'][index_result] \
                    = data['data'][field_variable_index]['components'][component_index]['values'][index_in]
            
            elif dimension == 3:
              for z in range(indices_begin[2],indices_end[2]):
                for y in range(indices_begin[1],indices_end[1]):
                  for x in range(indices_begin[0],indices_end[0]):
                    index_in = (z-indices_begin[2])*n_dofs_local[1]*n_dofs_local[0] \
                      + (y-indices_begin[1])*n_dofs_local[0] + (x-indices_begin[0])
                    index_result = z*n_dofs_global[1]*n_dofs_global[0] + y*n_dofs_global[0] + x
                    merged_data['data'][field_variable_index]['components'][component_index]['values'][index_result] \
                      = data['data'][field_variable_index]['components'][component_index]['values'][index_in]

    loaded_data.append(merged_data)
  return loaded_data


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
  
