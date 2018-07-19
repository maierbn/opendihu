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

# to make the loop easier
variation_euler=("N=3","N=4","N=5","N=6","N=7","N=8")
variation_heun=("N=3","N=4","N=5","N=6","N=7")
# the relative error for k=state_number(=57) states. One for euler, one
# for heun, and both contain information of different step sizes
err_rel_euler= {key : [None]*state_number for key in variation_euler}
err_rel_heun = {key : [None]*state_number for key in variation_heun}
err_abs_euler= {key : [None]*state_number for key in variation_euler}
err_abs_heun = {key : [None]*state_number for key in variation_heun}
data_euler={"N=3":data_euler_3,"N=4":data_euler_4,"N=5":data_euler_5,"N=6":data_euler_6,"N=7":data_euler_7,"N=8":data_euler_8}
data_heun ={"N=3":data_heun_3,"N=4":data_heun_4,"N=5":data_heun_5,"N=6":data_heun_6,"N=7":data_heun_7}
keyToIndex_e= {key : list(err_rel_euler.keys()).index(key) for key in variation_euler}
keyToIndex_h= {key : list(err_rel_heun.keys()).index(key) for key in variation_heun}
# width: N=x. height: e and h . depth: state
# Tensor tells which data points (error values) are usable (up = usable point). 0 means not usable, 1 means usable
Tensor_up_rel = [[[0 for x in range(len(variation_euler))] for y in range(2)] for z in range(state_number)] #anzusprechen mit "äußeres zuerst" 
Tensor_up_abs = [[[0 for x in range(len(variation_euler))] for y in range(2)] for z in range(state_number)]
exact_states=[]
                   
#def is_number(s):
#  try:
#    float(s)
#    return True
#  except ValueError:
#    return False

####################################################################################################################################
######### compute errors, track min and max values of them and set the markers for Tensor_up to indicate valid errors ##############
####################################################################################################################################
# euler relative error computation:
err_rel_min_e=1
err_rel_max_e=0
for key in variation_euler:
  data_temp=data_euler[key]
  for state in range(0,state_number):
    comp = (2*state) % (state_number-1)
    if state <= (state_number-1)/2:
      val=0
    else:
      val=1
    sum_e = 0
    for eval_point in range(0,evaluation_point_number):
      sol = solution[eval_point]['data'][1]['components'][comp]['values'][val] # this is u_i(t_n) where i is the state component and eval_point=n
      if math.isnan(sol):
        print("Warning: nan in solution! See state {}".format(state))
      if sol == 0:
        err_rel_euler[key][state]="DBZ"
        break
      appr = data_temp[eval_point]['data'][1]['components'][comp]['values'][val]
      if math.isnan(appr):
        err_rel_euler[key][state]="myNAN"
        break
      sum_e += math.pow((appr-sol)/sol,2)
      if sum_e > 1e30:
        err_rel_euler[key][state]="myNAN"
        #Tensor_up_rel[state][0][key]=0 # is already 0
        break
    err_rel_euler[key][state]=math.sqrt(sum_e/evaluation_point_number)
    Tensor_up_rel[state][0][keyToIndex_e[key]] = 1
    if err_rel_euler[key][state] > err_rel_max_e:
      err_rel_max_e=err_rel_euler[key][state]
    if err_rel_euler[key][state] < err_rel_min_e:
      if err_rel_euler[key][state]==0:
        if not state in exact_states:
          exact_states.append(state)
          print("Error for state {} is zero. Nothing to be done here.".format(state,key))
      else:
        err_rel_min_e=err_rel_euler[key][state]
    
####################################################################################################################################
# heun relative error computation:
err_rel_min_h=1
err_rel_max_h=0
for key in variation_heun:
  data_temp=data_heun[key]
  for state in range(0,state_number):
    comp = (2*state) % (state_number-1)
    if state <= (state_number-1)/2:
      val=0
    else:
      val=1
    sum_h = 0
    for eval_point in range(0,evaluation_point_number):
      sol = solution[eval_point]['data'][1]['components'][comp]['values'][val]
      if sol == 0:
        err_rel_heun[key][state]="DBZ"
        break
      appr = data_temp[eval_point]['data'][1]['components'][comp]['values'][val]
      if math.isnan(appr):
        err_rel_heun[key][state]="myNAN"
        break
      sum_h += math.pow((appr-sol)/sol,2)
      if sum_h > 1e100:
        err_rel_heun[key][state]="myNAN"
        break
    err_rel_heun[key][state]=math.sqrt(sum_h/evaluation_point_number)
    Tensor_up_rel[state][1][keyToIndex_h[key]] = 1
    if err_rel_euler[key][state] > err_rel_max_e:
      err_rel_max_e=err_rel_euler[key][state]
    if err_rel_euler[key][state] < err_rel_min_e:
      if not err_rel_euler[key][state]==0:
        err_rel_min_e=err_rel_euler[key][state]
        
####################################################################################################################################
err_rel_max=max(err_rel_max_e,err_rel_max_h)
err_rel_min=min(err_rel_min_e,err_rel_min_h)
####################################################################################################################################
# euler absolute error computation:
err_abs_min_e=1
err_abs_max_e=0
for key in variation_euler:
  data_temp=data_euler[key]
  for state in range(0,state_number):
    comp = (2*state) % (state_number-1)
    if state <= (state_number-1)/2:
      val=0
    else:
      val=1
    sum_e = 0
    for eval_point in range(0,evaluation_point_number):
      sol = solution[eval_point]['data'][1]['components'][comp]['values'][val]
      appr = data_temp[eval_point]['data'][1]['components'][comp]['values'][val]
      if math.isnan(appr):
        err_abs_euler[key][state]="myNAN"
        break
      sum_e += math.pow(appr-sol,2)
      if sum_e > 1e30:
        err_abs_euler[key][state]="myNAN"
        break
    err_abs_euler[key][state]=math.sqrt(sum_e/evaluation_point_number)
    Tensor_up_abs[state][0][keyToIndex_e[key]] = 1
    if err_abs_euler[key][state] > err_abs_max_e:
      err_abs_max_e=err_abs_euler[key][state]
    if err_abs_euler[key][state] < err_abs_min_e:
      if not err_abs_euler[key][state]==0:
        err_abs_min_e=err_abs_euler[key][state]

####################################################################################################################################
# heun absolute error computation:
err_abs_min_h=1
err_abs_max_h=0
for key in variation_heun:
  data_temp=data_heun[key]
  for state in range(0,state_number):
    comp = (2*state) % (state_number-1)
    if state <= (state_number-1)/2:
      val=0
    else:
      val=1
    sum_h = 0
    for eval_point in range(0,evaluation_point_number):
      sol = solution[eval_point]['data'][1]['components'][comp]['values'][val]
      appr = data_temp[eval_point]['data'][1]['components'][comp]['values'][val]
      if math.isnan(appr):
        err_abs_heun[key][state]="myNAN"
        break
      sum_h += math.pow(appr-sol,2)
      if sum_h > 1e30:
        err_abs_heun[key][state]="myNAN"
        break
    err_abs_heun[key][state]=math.sqrt(sum_h/evaluation_point_number)
    Tensor_up_abs[state][1][keyToIndex_h[key]] = 1
    if err_abs_heun[key][state] > err_abs_max_h:
      err_abs_max_h=err_abs_heun[key][state]
    if err_abs_heun[key][state] < err_abs_min_h:
      if not err_abs_heun[key][state]==0:
        err_abs_min_h=err_abs_heun[key][state]
        
####################################################################################################################################
err_abs_max=max(err_abs_max_e,err_abs_max_h)
err_abs_min=min(err_abs_min_e,err_abs_min_h)
####################################################################################################################################
####################################################################################################################################
#to plot the non-problematic states first, we need them sorted. If all erros for N=x exist and the convergence rate of Heun seems to
#be there, we put them into 'nice_states'. if the errors are there, but the rate is bad, we put them into 'problem_states'. Last, if
#the error doesnt exist for all N=x, we put them into 'late_convergence_states' (no 'late_problem_states' vector, yet).
#Note that there also exists the vector 'exact_states'.

def is_nice(vec):
  for x_j in vec:
    if x_j==0:
      return False
  return True

def is_converged(err_start,err_end,t_start,t_end):
  if err_end/err_start > 2*pow(t_end/t_start,2):
    return False
  else:
    return True
  
nice_states=[]
problem_states=[]
late_convergence_states=[]
for state in range(0,state_number):
  vec_to_test=[]
  if not state in exact_states:
    for key in variation_euler:
      vec_to_test.append(Tensor_up_abs[state][0][keyToIndex_e[key]])
    for key in variation_heun:
      vec_to_test.append(Tensor_up_abs[state][1][keyToIndex_h[key]])
    if is_nice(vec_to_test):
      if is_converged(err_abs_heun[variation_heun[0]][state],err_abs_heun[variation_heun[-1]][state],time_stepsizes_heun[0],time_stepsizes_heun[-1]):
        nice_states.append(state)
      else:
        print("No convergent behaviour in state {}.".format(state))
        problem_states.append(state)
    else:
      print("Late convergence might occur in state {}.".format(state))
      late_convergence_states.append(state)

# prepare plot
palette = plt.get_cmap('Set1')
#fig1 = plt.figure(1)
##plt.loglog(time_stepsizes_euler,[err_rel_euler[variation_euler[0]][0],err_rel_euler[variation_euler[1]][0],err_rel_euler[variation_euler[2]][0],err_rel_euler[variation_euler[3]][0]])
##margin = abs(max_value - min_value) * 0.1
#ax = plt.axes(xlim=(0.9*min(time_stepsizes_euler),1.1*max(time_stepsizes_euler)),ylim=(0.9*err_min,1.1*err_max))
#plt.xscale('log')#, basex=2)
#plt.grid(b=1,which='major',axis='both')
#plt.yscale('log')
#text = plt.figtext(0.15,0.16,"${err_i}^2 := \int_0^{100 ms}((u_i(t)-y_i(t))) / u_i(t) )^2 dt$",size=14)
#text_e = plt.figtext(0.15,0.83,"--x-- Euler",size=10)
#text_h = plt.figtext(0.15,0.78,"$-$+  Heun",size=10)
#title=plt.title("Relative errors depending on applied time step size")
#ax.set_xlabel('time step size')
#ax.set_ylabel('relative error')


## plot some lines:
#line, = ax.plot([], [], 'o-', lw=2)
#for ii in valid_rel_errors:
  #temp_e = [err_rel_euler[variation_euler[0]][ii],err_rel_euler[variation_euler[1]][ii],err_rel_euler[variation_euler[2]][ii],err_rel_euler[variation_euler[3]][ii],err_rel_euler[variation_euler[4]][ii],err_rel_euler[variation_euler[5]][ii]]
  #temp_h = [err_rel_heun[variation_heun[0]][ii],err_rel_heun[variation_heun[1]][ii],err_rel_heun[variation_heun[2]][ii],err_rel_heun[variation_heun[3]][ii],err_rel_heun[variation_heun[4]][ii]]
  #plt.plot(time_stepsizes_euler,temp_e,color=palette(ii % 7), marker='x', linestyle='--')
  #if temp_h[-1]/temp_h[0] > 2*pow(time_stepsizes_heun[-1]/time_stepsizes_heun[0],2):
    #plt.plot(time_stepsizes_heun,temp_h,color=palette(ii % 7), marker='o', linestyle='-')
  #else:
    #plt.plot(time_stepsizes_heun,temp_h,color=palette(ii % 7), marker='+', linestyle='-')
  #plt.plot([time_stepsizes_euler[-1], time_stepsizes_heun[-1]],[temp_e[-1], temp_h[-1]],color=palette(ii % 7), linestyle=':')#'k'
    
#plt.plot([max(time_stepsizes_euler),min(time_stepsizes_euler)],[err_max*0.8, err_max*0.8*min(time_stepsizes_euler)/max(time_stepsizes_euler)],'k-.')
#plt.plot([max(time_stepsizes_euler),min(time_stepsizes_euler)],[err_min*1.15/pow(min(time_stepsizes_euler)/max(time_stepsizes_euler),2), err_min*1.15],'k-.')
#print("Showing {} relative errors".format(len(valid_rel_errors)))

#plt.savefig(plot_name)

#if show_plot:
  #plt.show()


##ToDo: compute the absolute errors of the problematic states listed in the array 'problematic_states' and plot them to a second figure. They might provide more insight why there is no convergence.
#ps_number=len(problematic_states)
#if ps_number>0:
        
  #valid_abs_errors=[]
  #for ii in range(0,state_number):
    #if valid_abs_test(ii):
      #valid_abs_errors.append(ii)
    #else:
      #print("Problem in state {}, regarding absolute errors".format(ii))
  
  #fig2 = plt.figure(2)
  #ax = plt.axes(xlim=(0.9*min(time_stepsizes_euler),1.1*max(time_stepsizes_euler)),ylim=(0.9*err_abs_min,1.1*err_abs_max))
  #plt.xscale('log')#, basex=2)
  #plt.grid(b=1,which='major',axis='both')
  #plt.yscale('log')
  #text = plt.figtext(0.15,0.16,"${err_i}^2 := \int_0^{100 ms}(u_i(t)-y_i(t))^2 dt$",size=14)
  #text_e = plt.figtext(0.15,0.83,"--x-- Euler",size=10)
  #text_h = plt.figtext(0.15,0.78,"$-$+  Heun",size=10)
  #title=plt.title("Absolute errors depending on applied time step size")
  #ax.set_xlabel('time step size')
  #ax.set_ylabel('Error')

## plot some lines:
  #line, = ax.plot([], [], 'o-', lw=2)
  #for state in range(0,state_number):
    #temp_e=[]
    #temp_t_e=[]
    #temp_h=[]
    #temp_t_h=[]
    #M_usable=useable(state)
    #for err_j in range(0,len(M_usable[0])):
      #if M_usable[0][err_j]==1:
        #temp_e.append(err_abs_euler[variation_euler[err_j]][state])
        #temp_t_e.append(time_stepsizes_euler[err_j])
      #if M_usable[1][err_j]==1:
        #temp_h.append(err_abs_heun[variation_heun[err_j]][state])
        #temp_t_h.append(time_stepsizes_heun[err_j])
    #if valid_abs_test(state):
      #if not state in problematic_states:
        #plt.plot(temp_t_e,temp_e,color='grey', marker='x', linestyle='--')
        #plt.plot(temp_t_h,temp_h,color='grey', marker='+', linestyle='-')
        #plt.plot([temp_t_e[-1], temp_t_h[-1]],[temp_e[-1], temp_h[-1]],color='grey', linestyle=':')
      #else:
        #plt.plot(temp_t_e,temp_e,color=palette(state % 7), marker='x', linestyle='--')
        #plt.plot(temp_t_h,temp_h,color=palette(state % 7), marker='o', linestyle='-')
        #plt.plot([temp_t_e[-1], temp_t_h[-1]],[temp_e[-1], temp_h[-1]],color=palette(state % 7), linestyle=':')
    #else:
      ##print("time: {}".format(temp_t_h))
      ##print("error:{}".format(temp_h))
      #print("Convergence for state {} observed later.".format(state))
      #plt.plot(temp_t_e,temp_e,color=palette(state % 7), marker='x', linestyle='--', linewidth=1.0)
      #plt.plot(temp_t_h,temp_h,color=palette(state % 7), marker='+', linestyle='-', linewidth=1.0)
      #plt.plot([temp_t_e[-1], temp_t_h[-1]],[temp_e[-1], temp_h[-1]],color='grey', linestyle=':', linewidth=0.5)
        
      
  ##for ii in valid_abs_errors:
  ##  if not ii in problematic_states:
  ##    temp_e = [err_abs_euler[variation_euler[0]][ii],err_abs_euler[variation_euler[1]][ii],err_abs_euler[variation_euler[2]][ii],err_abs_euler[variation_euler[3]][ii],err_abs_euler[variation_euler[4]][ii],err_abs_euler[variation_euler[5]][ii]]
  ##    temp_h = [err_abs_heun[variation_heun[0]][ii],err_abs_heun[variation_heun[1]][ii],err_abs_heun[variation_heun[2]][ii],err_abs_heun[variation_heun[3]][ii],err_abs_heun[variation_heun[4]][ii]]
  ##    plt.plot(time_stepsizes_euler,temp_e,color='grey', marker='x', linestyle='--')
  ##    plt.plot(time_stepsizes_heun,temp_h,color='grey', marker='+', linestyle='-')
  ##    plt.plot([time_stepsizes_euler[-1], time_stepsizes_heun[-1]],[temp_e[-1], temp_h[-1]],color='grey', linestyle=':')#'k'
  ##for ii in problematic_states:
  ##  temp_e = [err_abs_euler[variation_euler[0]][ii],err_abs_euler[variation_euler[1]][ii],err_abs_euler[variation_euler[2]][ii],err_abs_euler[variation_euler[3]][ii],err_abs_euler[variation_euler[4]][ii],err_abs_euler[variation_euler[5]][ii]]
  ##  temp_h = [err_abs_heun[variation_heun[0]][ii],err_abs_heun[variation_heun[1]][ii],err_abs_heun[variation_heun[2]][ii],err_abs_heun[variation_heun[3]][ii],err_abs_heun[variation_heun[4]][ii]]
  ##  plt.plot(time_stepsizes_euler,temp_e,color=palette(ii % 7), marker='x', linestyle='--')
  ##  plt.plot(time_stepsizes_heun,temp_h,color=palette(ii % 7), marker='o', linestyle='-')
  ##  plt.plot([time_stepsizes_euler[-1], time_stepsizes_heun[-1]],[temp_e[-1], temp_h[-1]],color=palette(ii % 7), linestyle=':')#'k'
  ##  comp = (2*ii) % (state_number-1)
  ##  if ii <= (state_number-1)/2:
  ##    val=0
  ##  else:
  ##    val=1
  ##  print("Component {}: {}".format(ii,solution[evaluation_point_number-1]['data'][1]['components'][comp]['values'][val]))
    ##else:
    ##  temp_e = [err_abs_euler[variation_euler[0]][ii],err_abs_euler[variation_euler[1]][ii],err_abs_euler[variation_euler[2]][ii],err_abs_euler[variation_euler[3]][ii],err_abs_euler[variation_euler[4]][ii],err_abs_euler[variation_euler[5]][ii]]
    ##  temp_h = [err_abs_heun[variation_heun[0]][ii],err_abs_heun[variation_heun[1]][ii],err_abs_heun[variation_heun[2]][ii],err_abs_heun[variation_heun[3]][ii],err_abs_heun[variation_heun[4]][ii]]
    ##  plt.plot(time_stepsizes_euler,temp_e,color='grey', marker='x', linestyle='--')
    ##  plt.plot(time_stepsizes_heun,temp_h,color='grey', marker='+', linestyle='-')
    ##  plt.plot([time_stepsizes_euler[-1], time_stepsizes_heun[-1]],[temp_e[-1], temp_h[-1]],color=palette(ii % 7), linestyle=':')#'k'
  
  #plt.plot([max(time_stepsizes_euler),min(time_stepsizes_euler)],[err_abs_max*0.8, err_abs_max*0.8*min(time_stepsizes_euler)/max(time_stepsizes_euler)],'k-.')
  #plt.plot([max(time_stepsizes_euler),min(time_stepsizes_euler)],[err_abs_min*1.15/pow(min(time_stepsizes_euler)/max(time_stepsizes_euler),2), err_abs_min*1.15],'k-.')

  #plt.savefig("absoluteErrors.pdf")  
  #if show_plot:
    #plt.show()
  #print(err_abs_euler)
  #print(err_abs_heun)
  
#  for state in range(0,state_number):
#    print(useable(state))
print("Captured {} states in total (out of {}).".format(len(exact_states)+len(nice_states)+len(problem_states)+len(late_convergence_states),state_number))
print("  {} constant ones (thus, exact),".format(len(exact_states)))
print("  {} nicely converging ones,".format(len(nice_states)))
print("  {} problematic ones and".format(len(problem_states)))
print("  {} converging lately.".format(len(late_convergence_states)))
  
sys.exit(0)
