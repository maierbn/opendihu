#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# usage: ./merge_com_files.py <files>
#
# Merge the given com files to create a new file "output.com".

import sys

filenames = []
for filename in sys.argv[1:]:
  filenames.append(filename)
  
print filenames

current_time = 0

parsed_files = {}
exfile_base = ""

groups = set()

# parse files
for filename in filenames:
  
  files_started = False
  collected_lines = []
  last_lines_started = False
  last_lines = []
  invalid_dataset = False
  
  print "file \"{}\"".format(filename)
  with open(filename, 'rb') as f:
    for line in f.readlines():
      
      if "$fname =" in line:
        if files_started and not invalid_dataset:
            
          # add entry for current_time
          if current_time not in parsed_files:
            parsed_files[current_time] = []
          parsed_files[current_time].append(
            {"fname": exfile_base, "n_nodes": n_nodes, "n_elements": n_elements, "lines": collected_lines})
        collected_lines = []
          
        files_started = True
        invalid_dataset = False
        
        exfile_base = line[line.find("\"")+1:-2]
    
        # determine number of nodes
        try:
          with open(exfile_base+".exnode","rb") as exnode_f:
            for line2 in exnode_f.readlines():
              if "Node: " in line2:
                node_no = int(line2[7:])
          n_nodes = node_no
        except:
          invalid_dataset = True
          continue
        #print("n nodes: {}".format(n_nodes))
        
        # determine number of elements
        try:
          with open(exfile_base+".exelem","rb") as exelem_f:
            for line2 in exelem_f.readlines():
              if "Element: " in line2:
                substr = line2[21:]
                element_no = int(substr[:substr.find(" ")])
          n_elements = element_no
        except:
          invalid_dataset = True
          continue
        #print("n elements: {}".format(n_elements))
      
      if " time " in line:
        substr = line[line.find(" time ")+6:]
        current_time = float(substr[:substr.find(" ")])
      
      if "$group = " in line:
        group = line[line.find("$group = \"")+10:]
        group = group[:group.find("\"")]
        groups.add(group)
      
      if "$group = " in line or "# set the group name of the mesh" in line:
        last_lines_started = True      
        
      if last_lines_started:
        last_lines.append(line)
      else:
        collected_lines.append(line)
    if files_started:
        
      # add entry for current_time
      if current_time not in parsed_files:
        parsed_files[current_time] = []
      #print "collected_lines:",collected_lines
      parsed_files[current_time].append(
          {"fname": exfile_base, "n_nodes": n_nodes, "n_elements": n_elements, "lines": collected_lines})
        
# output everything
with open("output.com","wb") as outfile:
  outfile.write("# Merged com file by merge_com_files.py, \n# from files {}\n#\n#\n".format(filenames))
  for (time,d) in sorted(parsed_files.items()):
    outfile.write("# ------------------- time {} -------------------".format(time))
    
    node_offset = 1
    element_offset = 1
    
    for files_for_time in d:
      outfile.write("\n")
      #print "lines: [",files_for_time["lines"],"]"
      
      for line in files_for_time["lines"]:
        if "node_offset " in line:
          p0 = line.find("node_offset ")+12
          p1 = p0+line[p0:].find(" ")
          line = line[:p0] + str(node_offset) + line[p1:]
          
        if "element_offset " in line:
          p0 = line.find("element_offset ")+15
          p1 = p0+line[p0:].find(" ")
          line = line[:p0] + str(element_offset) + line[p1:]
          
        outfile.write(line)
      node_offset += files_for_time["n_nodes"]
      element_offset += files_for_time["n_elements"]
      
      
  for line in last_lines:
    if "gfx modify g_element $group points" in line:
      for group in groups:
        outfile.write("gfx modify g_element \""+group+"\""+line[27:])
    else:
      outfile.write(line)

  
  
  
