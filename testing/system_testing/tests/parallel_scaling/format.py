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
  
  t = c_longlong(int(value))
  
  first = True
  
  if t.value > c_longlong(3600000000000L).value:
    h = int(t.value / c_longlong(3600000000000L).value)    
    t.value -= h * c_longlong(3600000000000L).value
    
    if first:
      result += "{}h  ".format(h)
      first = False
    else:
      result += "{}h ".format(h)

  if t.value > c_longlong(60000000000L).value:
    minutes = int(t.value / c_longlong(60000000000L).value)
    t.value -= minutes * c_longlong(60000000000L).value
    
    if first:
      result += "{}min  ".format(minutes)
      first = False
    else:
      result += "{}min ".format(minutes)

  if t.value > c_longlong(1000000000L).value:    
    sec = int(t.value / c_longlong(1000000000L).value)
    t.value -= sec * c_longlong(1000000000L).value
    
    if first:
      result += "{}s  ".format(sec)
      first = False
    else:
      result += "{}s ".format(sec)

  if t.value > c_longlong(1000000L).value:
    msec = int(t.value / c_longlong(1000000L).value)
    t.value -= msec * c_longlong(1000000L).value
    
    if first:
      result += "{}ms  ".format(msec)
      first = False
    else:
      result += "{}ms ".format(msec)

  if t.value > 1000:
    usec = int(t.value / 1000)
    t.value -= usec * 1000
    
    if first:
      result += "{}us  ".format(usec)
      first = False
    else:
      result += "{}us ".format(usec)

  if (t > 0):
    
    if first:
      result += "{}ns  ".format(t.value)
      first = False
    else:
      result += "{}ns".format(t.value)

  #result += " (={}ns)".format(value);
  return result
    
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
      print "usage: {} <number> <format>\n Output number with a special format.\nValid values for <format> are:\n   mem   Format as number of bytes\n   ns   Format as duration in nanoseconds\n   n   Format as big number.".format(sys.argv[0])
      quit()
    
    number = int(sys.argv[1])
    
    formatSpec = "n"
    if len(sys.argv) == 3:
       formatSpec = sys.argv[2]
       
    if formatSpec == "mem":
      print str_format_memory(number)
    elif formatSpec == "ns":
      print str_format_nanoseconds(number)
    else:
      print str_format_number(number)
       
