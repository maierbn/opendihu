#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Transform C++ method signatures from header files to sphinx code that can be used for the documentation
#

import sys,os

if len(sys.argv) == 1:
  print("usage: {} <input file *.h>".format(sys.argv[0]))
  sys.exit(0)
  
input_filename = sys.argv[1]

debug = False

output = ""
def handle_method(method, comment):
  global output
  
  method_output = ""
  try:
    # parse method string
    pos = method.find("(")+1
    
    args = []
    end = False
    while not end:
      if "," in method[pos:]:
        next_pos = pos+method[pos:].find(",")
      elif ")" in method[pos:]:
        next_pos = pos+method[pos:].find(")")
        end = True
      else:
        next_pos = len(method)
        if debug: print("end")
        end = True
      
      if debug: print("part: {},{} [{}]".format(pos,next_pos,method[pos:next_pos]))
      
      new_arg = method[pos:next_pos].strip()
      
      added_arg = False
      if len(args) > 0:
        n_open = args[-1].count("<") - args[-1].count(">")
        if n_open > 0:
          args[-1] = args[-1]+","+new_arg
          added_arg = True
          
      if not added_arg:
        if new_arg != "":
          args.append(new_arg)
      pos = next_pos+1
      if debug: print("new pos: {}, end:{}".format(pos, end))

    if debug: print("")

    if ";" in method:
      method = method[0:method.find(";")]

    comment = comment.strip()
    if comment[0].islower():
      comment = comment[0].upper() + comment[1:]
    if comment[-1] != ".":
      comment += "."
    comment = comment.replace("_",'\_')
    
    if debug: 
      print("comment: [{}]".format(comment))
      print("args: {}".format(args))
      
    method_output += ".. cpp:function:: {}\n  \n  {}\n  \n".format(method, comment)

    if debug: 
      print("")
      print("")
      
    for arg in args:
      arg = arg.replace(":",'\:')
      
      method_output += "  :param {}: \n".format(arg)
    method_output += "  \n"
  except:
    if debug: 
      print("exception occured\n")
    
  output += method_output

# parse header file line-by-line
line_type = "none"
comment = ""
method = ""

with open(input_filename, "r") as f:
  lines = f.readlines()

for line in lines:
  line = line.strip()

  if debug: 
    print("line: [{}]".format(line))

  if "//!" in line:
    line_type = "comment"
    pos = line.find("//!")
    comment += line[pos+3:]
  
  elif len(line) == 0 or "}" in line:
    line_type = "empty"
    
    if ";" in method and "(" in method and ")" in method:
      handle_method(method, comment)
    
    # reset variables
    comment = ""
    method = ""
    
  else:
    line_type = "method"
    method += line

print("")
print(output)
