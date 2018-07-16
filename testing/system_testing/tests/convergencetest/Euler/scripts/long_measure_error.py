## # !/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Script to visualize python output files.
#

import sys
import numpy as np
import math

import csv
import collections
import copy
# from sets import Set # sets is already included in python3
import os
import time
import pickle
import py_reader    # reader utility for opendihu *.py files

files = ""

plot_name_condition = lambda filename: ".pdf" in filename

show_plot = True
if len(sys.argv) > 0:
  try:
    show_plot = int(sys.argv[1])
    files = sys.argv[2:]
  except:
    files = sys.argv[1:]
  plot_name_lst = list(np.extract(np.array(list(map(plot_name_condition, sys.argv))), sys.argv))
  if len(plot_name_lst)!=1:
    plot_name="long_fig.pdf"
  else:
    plot_name=plot_name_lst[0]
	
else:
  # get all input data in current directory
  ls = os.listdir(".")

  # sort files by number in file name
  files = sorted(ls)


# import needed packages from matplotlib
#if not show_plot:
#  import matplotlib as mpl
#  mpl.use('Agg')

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
from matplotlib import cm
#import matplotlib
from matplotlib.patches import Polygon

# right now, this is only implemented for a fix structure and a free but fix amount of input files.
# evaluation_point_number =  100
state_number = 57
# path information, fix
base_path="/usr/local/home/kraemer/software/opendihu/testing/system_testing/tests/convergencetest/"
#dir2="build_debug/out/N_2"
dir3="build_debug/out/N_3"
dir4="build_debug/out/N_4"
dir5="build_debug/out/N_5"
dir6="build_debug/out/N_6"
dir7="build_debug/out/N_7"
dir8="build_debug/out/N_8"
# read Euler Data of different stepsizes
#file_list_euler_2=[]
file_list_euler_3=[]
file_list_euler_4=[]
file_list_euler_5=[]
file_list_euler_6=[]
file_list_euler_7=[]
file_list_euler_8=[]
#for file in os.listdir(base_path+"Euler/"+dir2):
#  if file.endswith(".py"):
#    file_list_euler_2.append(base_path+"Euler/"+dir2+"/"+file)
#files = sorted(file_list_euler_2)
#data_euler_2=py_reader.load_data(files)
for file in os.listdir(base_path+"Euler/"+dir3):
  if file.endswith(".py"):
    file_list_euler_3.append(base_path+"Euler/"+dir3+"/"+file)    
if len(file_list_euler_3) == 0:
  print( "Data error. Did you run the tests before measurements?")
  sys.exit(0)
files = sorted(file_list_euler_3)
data_euler_3=py_reader.load_data(files)
for file in os.listdir(base_path+"Euler/"+dir4):
  if file.endswith(".py"):
    file_list_euler_4.append(base_path+"Euler/"+dir4+"/"+file)
files = sorted(file_list_euler_4)
data_euler_4=py_reader.load_data(files)
for file in os.listdir(base_path+"Euler/"+dir5):
  if file.endswith(".py"):
    file_list_euler_5.append(base_path+"Euler/"+dir5+"/"+file)
files = sorted(file_list_euler_5)
data_euler_5=py_reader.load_data(files)
for file in os.listdir(base_path+"Euler/"+dir6):
  if file.endswith(".py"):
    file_list_euler_6.append(base_path+"Euler/"+dir6+"/"+file)
files = sorted(file_list_euler_6)
data_euler_6=py_reader.load_data(files)
for file in os.listdir(base_path+"Euler/"+dir7):
  if file.endswith(".py"):
    file_list_euler_7.append(base_path+"Euler/"+dir7+"/"+file)
files = sorted(file_list_euler_7)
data_euler_7=py_reader.load_data(files)
for file in os.listdir(base_path+"Euler/"+dir8):
  if file.endswith(".py"):
    file_list_euler_8.append(base_path+"Euler/"+dir8+"/"+file)
files = sorted(file_list_euler_8)
data_euler_8=py_reader.load_data(files)
# read Heun Data of different stepsizes
#file_list_heun_2=[]
file_list_heun_3=[]
file_list_heun_4=[]
file_list_heun_5=[]
file_list_heun_6=[]
file_list_heun_7=[]
file_list_heun_8=[]
#for file in os.listdir(base_path+"Heun/"+dir2):
#  if file.endswith(".py"):
#    file_list_heun_2.append(base_path+"Heun/"+dir2+"/"+file)
#files = sorted(file_list_heun_2)
#data_heun_2=py_reader.load_data(files)
for file in os.listdir(base_path+"Heun/"+dir3):
  if file.endswith(".py"):
    file_list_heun_3.append(base_path+"Heun/"+dir3+"/"+file)
files = sorted(file_list_heun_3)
data_heun_3=py_reader.load_data(files)
for file in os.listdir(base_path+"Heun/"+dir4):
  if file.endswith(".py"):
    file_list_heun_4.append(base_path+"Heun/"+dir4+"/"+file)
files = sorted(file_list_heun_4)
data_heun_4=py_reader.load_data(files)
for file in os.listdir(base_path+"Heun/"+dir5):
  if file.endswith(".py"):
    file_list_heun_5.append(base_path+"Heun/"+dir5+"/"+file)
files = sorted(file_list_heun_5)
data_heun_5=py_reader.load_data(files)
for file in os.listdir(base_path+"Heun/"+dir6):
  if file.endswith(".py"):
    file_list_heun_6.append(base_path+"Heun/"+dir6+"/"+file)
files = sorted(file_list_heun_6)
data_heun_6=py_reader.load_data(files)
for file in os.listdir(base_path+"Heun/"+dir7):
  if file.endswith(".py"):
    file_list_heun_7.append(base_path+"Heun/"+dir7+"/"+file)
files = sorted(file_list_heun_7)
data_heun_7=py_reader.load_data(files)
for file in os.listdir(base_path+"Heun/"+dir8):
  if file.endswith(".py"):
    file_list_heun_8.append(base_path+"Heun/"+dir8+"/"+file)
files = sorted(file_list_heun_8)
data_heun_8=py_reader.load_data(files)
# use highest resolution Heun method as solution
solution=data_heun_8
# extrakt time step information
time_stepsizes_euler=[data_euler_3[-1]["currentTime"]/data_euler_3[-1]["timeStepNo"],data_euler_4[-1]["currentTime"]/data_euler_4[-1]["timeStepNo"],data_euler_5[-1]["currentTime"]/data_euler_5[-1]["timeStepNo"],data_euler_6[-1]["currentTime"]/data_euler_6[-1]["timeStepNo"],data_euler_7[-1]["currentTime"]/data_euler_7[-1]["timeStepNo"],data_euler_8[-1]["currentTime"]/data_euler_8[-1]["timeStepNo"]]
time_stepsizes_heun=[data_heun_3[-1]["currentTime"]/data_heun_3[-1]["timeStepNo"],data_heun_4[-1]["currentTime"]/data_heun_4[-1]["timeStepNo"],data_heun_5[-1]["currentTime"]/data_heun_5[-1]["timeStepNo"],data_heun_6[-1]["currentTime"]/data_heun_6[-1]["timeStepNo"],data_heun_7[-1]["currentTime"]/data_heun_7[-1]["timeStepNo"]]
time_stepsize_solution=[data_heun_8[-1]["currentTime"]/data_heun_8[-1]["timeStepNo"]]

# get number evaluation points (output files per setting)
# must be the same for all settings!
evaluation_point_number=len(file_list_euler_3)
epn=evaluation_point_number
if (len(file_list_euler_3)!=epn or len(file_list_euler_4)!=epn or len(file_list_euler_5)!=epn or len(file_list_euler_6)!=epn or len(file_list_euler_7)!=epn or len(file_list_euler_8)!=epn or len(file_list_heun_3)!=epn or len(file_list_heun_4)!=epn or len(file_list_heun_5)!=epn or len(file_list_heun_6)!=epn or len(file_list_heun_7)!=epn):
  print("WARNING: Input file numbers differ for different time step sizes!")
  print("         >> The computation of the relative error might be wrong.")
  print("")

# the relative error for k=state_number=57 states. One for euler, one
# for heun, and both contain information of different step sizes
err_rel_euler = {"N=3" : [None] * state_number,"N=4" : [None] * state_number,"N=5" : [None] * state_number,"N=6" : [None] * state_number,"N=7" : [None] * state_number,"N=8" : [None] * state_number}
err_rel_heun = {"N=3" : [None] * state_number,"N=4" : [None] * state_number, "N=5" : [None] * state_number,"N=6" : [None] * state_number,"N=7" : [None] * state_number}

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def valid_test(index):
  if (is_number(err_rel_euler["N=3"][index]) and is_number(err_rel_euler["N=4"][index]) and is_number(err_rel_euler["N=5"][index])and is_number(err_rel_euler["N=6"][index]) and is_number(err_rel_euler["N=7"][index]) and is_number(err_rel_euler["N=8"][index])) :
    if (is_number(err_rel_heun["N=3"][index]) and is_number(err_rel_heun["N=4"][index])and is_number(err_rel_heun["N=5"][index]) and is_number(err_rel_heun["N=6"][index]) and is_number(err_rel_heun["N=7"][index])) :
      return True
  else:
    return False
    
# to make the loop easier
data_euler={"N=3":data_euler_3,"N=4":data_euler_4,"N=5":data_euler_5,"N=6":data_euler_6,"N=7":data_euler_7,"N=8":data_euler_8}
data_heun={"N=3":data_heun_3,"N=4":data_heun_4,"N=5":data_heun_5,"N=6":data_heun_6,"N=7":data_heun_7}
variation_euler=("N=3","N=4","N=5","N=6","N=7","N=8")
# euler relative error computation:
for key in variation_euler:
  data_temp=data_euler[key]
  for state in range(0,state_number):
    sum_e = 0
    for eval_point in range(0,evaluation_point_number):
      sol = solution[eval_point]['data'][1]['components'][state]['values'][0]
      if math.isnan(sol):
        print("Warning: nan in solution! See state {}".format(state))
      if sol == 0:
        err_rel_euler[key][state]="DBZ"
        break
      appr = data_temp[eval_point]['data'][1]['components'][state]['values'][0]
      if math.isnan(appr):
        err_rel_euler[key][state]="myNAN"
        break
      sum_e += math.pow((appr-sol)/sol,2)
      if sum_e > 1e30:
        err_rel_euler[key][state]="myNAN"
        break
    else:
      err_rel_euler[key][state]=math.sqrt(sum_e/evaluation_point_number)
    
variation_heun=("N=3","N=4","N=5","N=6","N=7")
# heun relative error computation:
for key in variation_heun:
  data_temp=data_heun[key]
  for state in range(0,state_number):
    #if state==56:
    #  print("state=56")
    sum_h = 0
    for eval_point in range(0,evaluation_point_number):
      sol = solution[eval_point]['data'][1]['components'][state]['values'][0]
      if math.isnan(sol):
        print("Warning: nan in solution! See state {}".format(state))
    #  if state==56:
    #    print("sol={}".format(sol))
      if sol == 0:
        err_rel_heun[key][state]="DBZ"
        break
      appr = data_temp[eval_point]['data'][1]['components'][state]['values'][0]
      if math.isnan(appr):
        err_rel_heun[key][state]="myNAN"
        break
    #  if state==56:
    #    print("appr={}".format(appr))
    #  if eval_point==evaluation_point_number-1:
      sum_h += math.pow((appr-sol)/sol,2)
    #  if state==56:
    #    print("sum={}".format(sum_h))
      if sum_h > 1e100:
        err_rel_heun[key][state]="myNAN"
        break
    else:
      err_rel_heun[key][state]=math.sqrt(sum_h/evaluation_point_number)
    #  if state==56:
    #    print("err_rel={}".format(err_rel_heun[key][state]))
        
#print(err_rel_euler)
#print(err_rel_heun)

# skip all 0 and NAN errors if one of err_rel_* has one of those.
#I.e.: build an array [0,1,2,5,6,9,...,52] where e.g. 3,4,7,8,...,51,53,... are excluded
valid_errors=[]
for ii in range(0,state_number):
  if valid_test(ii):
    valid_errors.append(ii)
  else:
    print("Problem in state {}.".format(ii))

#min_value = min(dataVm)
#max_value = max(dataVm)
#min_x = time[0]
#max_x = time[-1]

#print( "value range: [{}, {}]".format(min_value, max_value))

# prepare plot
palette = plt.get_cmap('Set1')

err_min=1
err_max=0
for ii in valid_errors:
  if err_rel_heun[variation_heun[4]][ii]<err_min:
    err_min=err_rel_heun[variation_heun[4]][ii]
    for jj in range(0,4): # to capture oscillations in small errors
      if err_rel_heun[variation_heun[jj]][ii]<err_min:
        err_min=err_rel_heun[variation_heun[jj]][ii]
  if err_rel_euler[variation_euler[0]][ii]>err_max:
    err_max=err_rel_euler[variation_euler[0]][ii]

fig1 = plt.figure(1)
#plt.loglog(time_stepsizes_euler,[err_rel_euler[variation_euler[0]][0],err_rel_euler[variation_euler[1]][0],err_rel_euler[variation_euler[2]][0],err_rel_euler[variation_euler[3]][0]])
#margin = abs(max_value - min_value) * 0.1
ax = plt.axes(xlim=(0.9*min(time_stepsizes_euler),1.1*max(time_stepsizes_euler)),ylim=(0.9*err_min,1.1*err_max))
plt.xscale('log')#, basex=2)
plt.grid(b=1,which='major',axis='both')
plt.yscale('log')
text = plt.figtext(0.15,0.16,"${err_i}^2 := \int_0^{100 ms}((u_i(t)-y_i(t))) / u_i(t) )^2 dt$",size=14)
text_e = plt.figtext(0.15,0.78,"--x-- Euler",size=10)
text_h = plt.figtext(0.15,0.73,"$-$+  Heun",size=10)
title=plt.title("Relative errors depending on applied time step size")
ax.set_xlabel('time step size')
ax.set_ylabel('relative error')

problematic_states=[]

# plot some lines:
line, = ax.plot([], [], 'o-', lw=2)
for ii in valid_errors:
  temp_e = [err_rel_euler[variation_euler[0]][ii],err_rel_euler[variation_euler[1]][ii],err_rel_euler[variation_euler[2]][ii],err_rel_euler[variation_euler[3]][ii],err_rel_euler[variation_euler[4]][ii],err_rel_euler[variation_euler[5]][ii]]
  temp_h = [err_rel_heun[variation_heun[0]][ii],err_rel_heun[variation_heun[1]][ii],err_rel_heun[variation_heun[2]][ii],err_rel_heun[variation_heun[3]][ii],err_rel_heun[variation_heun[4]][ii]]
  plt.plot(time_stepsizes_euler,temp_e,color=palette(ii % 7), marker='x', linestyle='--')
  if temp_h[-1]/temp_h[0] > 2*pow(time_stepsizes_heun[-1]/time_stepsizes_heun[0],2):
    print("No convergent behaviour in state {}.".format(ii))
    problematic_states.append(ii)
    plt.plot(time_stepsizes_heun,temp_h,color=palette(ii % 7), marker='o', linestyle='-')
  else:
    plt.plot(time_stepsizes_heun,temp_h,color=palette(ii % 7), marker='+', linestyle='-')
  plt.plot([time_stepsizes_euler[-1], time_stepsizes_heun[-1]],[temp_e[-1], temp_h[-1]],color=palette(ii % 7), linestyle=':')#'k'
    
plt.plot([max(time_stepsizes_euler),min(time_stepsizes_euler)],[err_max*0.8, err_max*0.8*min(time_stepsizes_euler)/max(time_stepsizes_euler)],'k-.')
plt.plot([max(time_stepsizes_euler),min(time_stepsizes_euler)],[err_min*1.15/pow(min(time_stepsizes_euler)/max(time_stepsizes_euler),2), err_min*1.15],'k-.')
print("Showing {} relative errors".format(len(valid_errors)))

plt.savefig(plot_name)

if show_plot:
  plt.show()

#ToDo: compute the absolute errors of the problematic states listed in the array 'problematic_states' and plot them to a second figure. They might provide more insight why there is no convergence.
ps_number=len(problematic_states)
if ps_number>0:
  err_abs_euler = {"N=3" : [None] * state_number,"N=4" : [None] * state_number,"N=5" : [None] * state_number,"N=6" : [None] * state_number,"N=7" : [None] * state_number,"N=8" : [None] * state_number}
  err_abs_heun = {"N=3" : [None] * state_number,"N=4" : [None] * state_number, "N=5" : [None] * state_number,"N=6" : [None] * state_number,"N=7" : [None] * state_number}

  for key in variation_euler:
    data_temp=data_euler[key]
    for state in valid_errors:
      sum_e = 0
      for eval_point in range(0,evaluation_point_number):
        sol = solution[eval_point]['data'][1]['components'][state]['values'][0]
        appr = data_temp[eval_point]['data'][1]['components'][state]['values'][0]
        sum_e += math.pow((appr-sol),2)
        if sum_e > 1e90:
          err_abs_euler[key][state]="myNAN"
          break
      else:
        err_abs_euler[key][state]=math.sqrt(sum_e/evaluation_point_number)
    
  variation_heun=("N=3","N=4","N=5","N=6","N=7")
# heun relative error computation:
  for key in variation_heun:
    data_temp=data_heun[key]
    for state in valid_errors:
      sum_h = 0
      for eval_point in range(0,evaluation_point_number):
        sol = solution[eval_point]['data'][1]['components'][state]['values'][0]
        appr = data_temp[eval_point]['data'][1]['components'][state]['values'][0]
        #print(appr)
        sum_h += math.pow((appr-sol),2)
        #print(sum_h)
        if sum_h > 1e100:
          err_abs_heun[key][state]="myNAN"
          #print("break")
          break
      err_abs_heun[key][state]=math.sqrt(sum_h/evaluation_point_number)
      #print("added {}".format(math.sqrt(sum_h/evaluation_point_number)))
   
  err_min=1
  err_max=0
  for ii in problematic_states:
    if err_abs_heun[variation_heun[4]][ii]<err_min:
      err_min=err_abs_heun[variation_heun[4]][ii]
      for jj in range(0,4): # to capture oscillations in small errors
        if err_abs_heun[variation_heun[jj]][ii]<err_min:
          err_min=err_abs_heun[variation_heun[jj]][ii]
    if err_abs_euler[variation_euler[0]][ii]>err_max:
      err_max=err_abs_euler[variation_euler[0]][ii]
  
  fig2 = plt.figure(2)
  ax = plt.axes(xlim=(0.9*min(time_stepsizes_euler),1.1*max(time_stepsizes_euler)),ylim=(0.9*err_min,1.1*err_max))
  plt.xscale('log')#, basex=2)
  plt.grid(b=1,which='major',axis='both')
  plt.yscale('log')
  text = plt.figtext(0.15,0.16,"${err_i}^2 := \int_0^{100 ms}((u_i(t)-y_i(t))) / u_i(t) )^2 dt$",size=14)
  text_e = plt.figtext(0.15,0.78,"--x-- Euler",size=10)
  text_h = plt.figtext(0.15,0.73,"$-$+  Heun",size=10)
  title=plt.title("Absolute errors depending on applied time step size")
  ax.set_xlabel('time step size')
  ax.set_ylabel('Error')

# plot some lines:
  line, = ax.plot([], [], 'o-', lw=2)
  for ii in valid_errors:
    if ii in problematic_states:
      temp_e = [err_abs_euler[variation_euler[0]][ii],err_abs_euler[variation_euler[1]][ii],err_abs_euler[variation_euler[2]][ii],err_abs_euler[variation_euler[3]][ii],err_abs_euler[variation_euler[4]][ii],err_abs_euler[variation_euler[5]][ii]]
      temp_h = [err_abs_heun[variation_heun[0]][ii],err_abs_heun[variation_heun[1]][ii],err_abs_heun[variation_heun[2]][ii],err_abs_heun[variation_heun[3]][ii],err_abs_heun[variation_heun[4]][ii]]
      plt.plot(time_stepsizes_euler,temp_e,color=palette(ii % 7), marker='x', linestyle='--')
      plt.plot(time_stepsizes_heun,temp_h,color=palette(ii % 7), marker='o', linestyle='-')
      plt.plot([time_stepsizes_euler[-1], time_stepsizes_heun[-1]],[temp_e[-1], temp_h[-1]],color=palette(ii % 7), linestyle=':')#'k'
    else:
      temp_e = [err_abs_euler[variation_euler[0]][ii],err_abs_euler[variation_euler[1]][ii],err_abs_euler[variation_euler[2]][ii],err_abs_euler[variation_euler[3]][ii],err_abs_euler[variation_euler[4]][ii],err_abs_euler[variation_euler[5]][ii]]
      temp_h = [err_abs_heun[variation_heun[0]][ii],err_abs_heun[variation_heun[1]][ii],err_abs_heun[variation_heun[2]][ii],err_abs_heun[variation_heun[3]][ii],err_abs_heun[variation_heun[4]][ii]]
      plt.plot(time_stepsizes_euler,temp_e,color='grey', marker='x', linestyle='--')
      plt.plot(time_stepsizes_heun,temp_h,color='grey', marker='+', linestyle='-')
      plt.plot([time_stepsizes_euler[-1], time_stepsizes_heun[-1]],[temp_e[-1], temp_h[-1]],color=palette(ii % 7), linestyle=':')#'k'
  
  plt.plot([max(time_stepsizes_euler),min(time_stepsizes_euler)],[err_max*0.8, err_max*0.8*min(time_stepsizes_euler)/max(time_stepsizes_euler)],'k-.')
  plt.plot([max(time_stepsizes_euler),min(time_stepsizes_euler)],[err_min*1.15/pow(min(time_stepsizes_euler)/max(time_stepsizes_euler),2), err_min*1.15],'k-.')

  plt.savefig("plot_name.fig")  
  if show_plot:
    plt.show()
  #print(err_abs_euler)
  #print(err_abs_heun)
   
sys.exit(0)
