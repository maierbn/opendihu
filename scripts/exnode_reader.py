#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy as np
import matplotlib.pyplot as plt
import csv
import collections
import copy
from sets import Set
import os
import collections

geometry_list = []
vm_list = []

"""
  Parse .exnode file and extract the values for specified fields for each node.
  fields should be a list of items ["field_name", component_no], e.g. [["Vm",1],["Coordinate",1],["Coordinate",2]], or [["Vm",1]] (note the double squared brackets)
  The field specifier must exists in the exnode file. 
  Returns a dict: node=parse_file_dict(filename, fields), where node[node_no] is a list of the requested fields for node with that node number.
  Note that node_no starts with 1, as in the exnode file.
  You should not use this function, but parse_file, which returns a numpy array and also considers files from other processes.
"""
def parse_file_dict(filename, fields):
  
  with  open(filename) as f:
    
    #print "file ({})".format(filename)
    
    index_current_node = None
    next_line_contains_index = [False for _ in range(len(fields))]
    component_no = None
    line_indices = [-1 for _ in range(len(fields))]
    field_value = np.zeros((len(fields)))
    field_no = 0
    node = {}
    node_no = None
    
    for topline in f:
      
      # split lines with multiple values to lines that only contain one value
      sublines = [topline]
      if "e+" in topline or "e-" in topline:   # if the line contains values, e.g. "1.785880e-01"
        sublines = topline.split(" ")
      for line in sublines:
        if line.strip() == "":
          continue
          
        if component_no is not None:
          component_no += 1
        if index_current_node is not None:
          index_current_node += 1
        
        #print "line=[{}], component_no: {}, index_current_node: {}".format(line, component_no, index_current_node)
        
        # parse line index
        if component_no is not None:
          for i in range(len(fields)):
            if next_line_contains_index[i] and component_no == fields[i][1]:
              pos = line.index("index=")
              posend = line.index(",", pos)
              str = line[pos+len("index="):posend].strip()
              #print("{}:{}: [{}]".format(pos,posend,str))
              line_index = int(str)
              next_line_contains_index[i] = False
              #print "found index for \"{}\" component {}: {} (field_no={})".format(fields[i][0], fields[i][1], line_index, i)
              line_indices[i] = line_index
        
        # parse value of a field if the respective line is the current
        # loop over every field that should be extracted
        for i in range(len(fields)):
          if index_current_node == line_indices[i]:   # if the index of the values in the current node equals the index for the field to be extracted
            value = float(line.strip())     # remove whitespace and convert value to float
            field_value[i] = value          # store value
            #print("parse {line}, set \"{field}\"[{component}] = {value}, line index: {index}".format(line=line.strip(), field=fields[i][0], component=fields[i][1], value=value, index=line_indices[i]))
            
        # determine if the next line contains the statement of an index of a field (which happens in the header of the exnode file)
        for (i,field) in enumerate(fields):
          if ") "+field[0] in line:     # if ") fieldname" appears in the line, the next line(s) contain indices
            field_no = i
            next_line_contains_index[i] = True    # store the information that the next line contains the respective index statement
            component_no = 0                      # reset the counter for components for more index statements of the components of the field
            #print("field [{}] starts".format(field))
        
        if "Node:" in line:                 # if a new node begins in the current line
          component_no = None
          
          # save values of previous node
          if node_no != None:               # if there was a previous node
                          
            # assign the previously extracted values to the node
            node[node_no] = field_value
          
            # reset field_value for next node
            field_value = [0 for _ in range(len(field_value))]
            
          # parse node number of current node
          pos = line.find("Node:")+len("Node:")
          node_no = int(line[pos:].strip())
          #print("node no {} starts".format(node_no))
                    
          # reset the line counter of the current node
          index_current_node = 0
      
  # save value of last node
  if node_no != None:               # if there was a previous node
              
    # assign the previously extracted values to the node
    node[node_no] = field_value
  
    # reset field_value for next node
    field_value = [0 for _ in range(len(field_value))]
  
  return node
     
"""
  Parse .exnode file and extract the values for specified fields for each node.
  If a file like "Time_2_123.part0.exnode" is given also all other nodes with part1, part2 etc. are parsed automatically.
  "fields" should be a list of items ["field_name", component_no], e.g. [["Vm",1],["Coordinate",1],["Coordinate",2]], or [["Vm",1]] (note the double squared brackets)
  The field specifier must exists in the exnode file. 
  Return value is a numpy array where each row is one node and each column corresponds to one of the specified fields.
  Example: node=parse_file(filename, fields), then node[node_no-1,:] contains the values of the requested fields for the node with that node number.
  Note that the array index starts with 0 but the node indices start with 1.
  
  If a node between 1 and the maximum node number does not appear in the given file or its associated files from the other processes, zeros are inserted for that missing node.
"""
def parse_file(filename, fields):
  if not os.path.exists(filename):
    return None
  
  if ".part" in filename:
    pos = filename.find(".part")
    prefix = filename[0:pos]
    suffix = filename[filename.find(".",pos+5):]
    #print "filename[{}], prefix[{}], suffix[{}]".format(filename,prefix,suffix)
    
    # check how many files from different processes there are
    n_files = 0
    while True:
      if not os.path.isfile(prefix+".part"+str(n_files)+suffix):
        break
      n_files += 1
      
    # loop over files and accumulate the result
    result = {}
    for i in range(n_files):
      filename = prefix+".part"+str(i)+suffix
      data_single_file = parse_file_dict(filename, fields)
      
      result.update(data_single_file)
    
    # insert dummy values for nodes that were not present in the file(s)
    dummy = [0.0] * len(result.values()[0])
    number_of_nodes = max(result)
    missing_keys = set(range(1, number_of_nodes+1)).difference(result)
    if len(missing_keys) != 0:
      print "Warning: files \"{}\" do not contains values for {} node(s), beginning with node {}".\
      format(prefix+".part*"+suffix, len(missing_keys), list(missing_keys)[0])
    result.update((key, dummy) for key in missing_keys)
    
  else:
    result = parse_file_dict(filename, fields)
    
  # convert dictionary to numpy array
  ordered_result = collections.OrderedDict(result)
  sorted_result = sorted(ordered_result.items())
  result = np.array([item[1] for item in sorted_result])
    
  return result
  
    
    
"""
  Compare function for sort
  Compare OpenCMISS output files by their iteration number *_<no>_.part*
  This is not lexicographic ordering because of e.g. _1_, _2_, _10_, and not _1_, _10_, _2_
  This is a helper function for sort and should not be used directly.
"""
def compare(x,y):
  if ".part" in x:
    posend = x.find(".part")
    pos = x.rfind("_", 0, posend)+1
    strx = x[pos:posend]

  if ".part" in y:
    posend = y.find(".part")
    pos = y.rfind("_", 0, posend)+1
    stry = y[pos:posend]
  
  try:  
    nx = int(strx)
  except:
    nx = -5
    
  try:
    ny = int(stry)
  except:
    ny = -5
    
  return nx - ny
  
"""
  Get all the *.exnode files in the given directory and return the filenames as a list sorted by timestep in the filename
"""
def get_exnode_files(directory):
  # get all input data in current directory
  ls = os.listdir(directory)
  
  # sort files by number in file name
  directory_content = sorted(ls, cmp=compare)
  
  # extract the files that are exnode files
  exnode_condition = lambda filename: ".exnode" in filename
  directory_content = list(np.extract(map(exnode_condition, directory_content), directory_content))

  return directory_content

def parse_exelem_file(filename):

  with open(filename) as f:    
    #print "file ({})".format(filename)
    
    element_no = None
    next_line_contains_nodes = False
    element_nodes = {}
    
    for line in f:
      if "Element:" in line:
        pos = line.index("Element:")
        str = line[pos+len("Element:"):-1].strip()
        posend = str.index(" ")
        str = str[0:posend].strip()
        element_no = int(str)
        if element_no == 0:
          continue
        
        if element_no not in element_nodes:
          element_nodes[element_no] = []
        
      if next_line_contains_nodes:
        next_line_contains_nodes = False
        # parse node list
        for node_no_str in line.strip().split(" "):
          try:
            node_no = int(node_no_str.strip())
            element_nodes[element_no].append(node_no)
          except:
            pass
        
      if "Nodes:" in line:
        next_line_contains_nodes = True
    
  return element_nodes
    
if __name__ == "__main__":
  # execute only if run as a script
  
  filename_exnode = "Time_2_100.part0.exnode"
  if len(sys.argv) >= 2:
    filename_exnode = sys.argv[1]
    
    filename_exelem,_ = os.path.splitext(filename_exnode)
    filename_exelem += ".exelem"
    
    if len(sys.argv) >= 3:
      filename_exelem = sys.argv[2]
  else:
    print("usage: ./exnode_reader.py <filename.exnode> [<filename.exelem>]")
    sys.exit(0)
  
  print("filename exnode: \"{}\"".format(filename_exnode))
  print("filename exelem: \"{}\"".format(filename_exelem))
  
  nodal_values = parse_file(filename_exnode, [["coordinates",1],["coordinates",2],["coordinates",3]])
  
  print("result:")
  print(nodal_values)

  element_nodes = parse_exelem_file(filename_exelem)
  print("element_nodes:")
  print(element_nodes)
  
  
