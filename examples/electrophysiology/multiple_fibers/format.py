#!/usr/bin/env python

from ctypes import c_longlong
import math
import sys

def str_format_memory(value):
  '''
  format a value representing a number of bytes
  '''
  result = ""
  if value >= (1024 * 1024 * 1024 * 1024):
    tera = value / (1024. * 1024 * 1024 * 1024)
    result += "{:.3f} TiB".format(int(tera))

  elif (value >= (1024 * 1024 * 1024)):
    giga = value / (1024. * 1024 * 1024)
    result += "{:.3f} GiB".format(int(giga))

  elif (value >= (1024 * 1024)):
    mega = value / (1024. * 1024)
    result += "{:.3f} MiB".format(int(mega))

  elif (value >= (1024)):
    kilo = value / (1024.);
    result += "{:.3f} kiB".format(int(kilo))

  else:
    result += "{} B".format(int(value))

  if value >= 1024:
    result += " ("+str_format_number(int(value))+" Bytes)";
  
  return result
    
def str_format_nanoseconds(value):
  '''
  format a value representing a number of nanoseconds
  '''
  result = ""
  
  return value  
    
def str_format_seconds(value):
  '''
  format a float value representing a number of seconds
  '''
  #return str_format_nanoseconds(value*1000000000.0)
  return "{:0.5}".format(float(value))+" s"
    
def str_format_number(value):
  '''
  format a value representing a large number
  '''
  result = ""
  number = "{}".format(value)

  pos = len(number) % 3
  while pos < len(number):

    start_pos = max(0, (int)(pos - 3))
    
    if pos > 0:
      result += number[start_pos:pos]+"'";
    
    pos += 3

  result += number[max(0, (int)(pos - 3)):]
  return result;
    
if __name__ == '__main__':
    
    # read in arguments
    if len(sys.argv) == 1:
      print("usage: {} <number> <format>\n Output number with a special format.\nValid values for <format> are:\n   mem   Format as number of bytes\n   ns   Format as duration in nanoseconds\n   n   Format as big number.".format(sys.argv[0]))
      quit()
    
    number = int(sys.argv[1])
    
    formatSpec = "n"
    if len(sys.argv) == 3:
       formatSpec = sys.argv[2]
       
    if formatSpec == "mem":
      print(str_format_memory(number))
    elif formatSpec == "ns":
      print(str_format_nanoseconds(number))
    else:
      print(str_format_number(number))
       
