#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script to visualize python output files.
#
# usage: 
#   plot                                            # plot all python files in current directory
#   plot file1.py file2.py ...                      # plot specific files
#   plot monodomain*                                # plot all files starting with monodomain
#   plot 0 [<files.py>]                             # don't show plot window
#   plot membrane/V_s membrane/V_d  [<files.py>]    # plot the specified fields

import sys
import numpy as np

import csv
import collections
import copy
# from sets import Set # sets included in python3
import os
import time
import pickle
import py_reader    # reader utility for opendihu *.py files

filenames = ""
specified_field_names = []

show_field_variable_names = True
show_plot = True
if len(sys.argv) > 1:
  try:
    show_plot = int(sys.argv[1])
    filenames = sys.argv[2:]
  except:
    filenames = sys.argv[1:]
    
  for filename in filenames:
    if ".py" in filename:
      break
    else:
      specified_field_names.append(filename)
else:
  # get all input data in current directory
  ls = os.listdir(".")

  # sort file_names by number in file name
  filenames = sorted(ls)


# import needed packages from matplotlib
import matplotlib as mpl
if not show_plot:
  mpl.use('Agg')

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
from matplotlib import cm
from matplotlib.patches import Polygon
import matplotlib.gridspec as gridspec

# extract the filenames that are python files
solution_py_condition = lambda filename: ".py" in filename
solution_py_files = list(np.extract(np.array(list(map(solution_py_condition, filenames))), filenames))


# sort files by number in file name
solution_py_files = sorted(solution_py_files)

# output file names
if len(solution_py_files) == 1:
  print("1 file: {}".format(solution_py_files[0]))
elif len(solution_py_files) <= 4:
  print("{} files: {}".format(len(solution_py_files),solution_py_files))
else:
  print("{} files: {}, {}, {}, ..., {}".format(len(solution_py_files), solution_py_files[0],solution_py_files[1],solution_py_files[2],solution_py_files[-1]))

# load data
data = py_reader.load_data(solution_py_files)

if show_field_variable_names and data:  
  print("")
  print("Available field variables (use `plot <field variable or component names> <files>` to plot selected variables over time):")
  field_variable_names = py_reader.get_field_variable_names(data[0])
  for field_variable_name in field_variable_names:
    component_names = py_reader.get_component_names(data[0], field_variable_name)
    print("*  {}, {} components: {}".format(field_variable_name, len(component_names), component_names))
  print("")

if len(data) == 0:
  print( "no data found.")
  sys.exit(0)

# set global parameters for font sizes
plt.rcParams.update({'font.size': 14})
plt.rcParams['lines.linewidth'] = 3
plt.rcParams['lines.markersize'] = 8

dimension = data[0]['dimension']

####################
# 1D
if dimension == 1:
  
  fig = plt.figure(figsize=(10,8))
  
  show_geometry = True     # if the fibre geometry should be displayed in a 3D plot in a separate axis (ax2) on top of the solution plot
  show_components = False  # if all the components of the solution should be displayed
  show_specified_fields = False   # if only the specified fields should be plotted in 2D in a single axis
  plot_over_time = data[0]['nElements'] == [0]   # if the plot should have time as x-axis instead of geometry
  compute_stress = False
  
  if len(specified_field_names) > 0:
    show_specified_fields = True
    plot_over_time = True
  
  min_x, max_x = py_reader.get_min_max(data, "geometry", "x")
  min_y, max_y = py_reader.get_min_max(data, "geometry", "y")
  min_z, max_z = py_reader.get_min_max(data, "geometry", "z")
  if abs(min_x-max_x) < 1e-5 and abs(min_y-max_y) < 1e-5 and abs(min_z-max_z) < 1e-5:
    plot_over_time = True

  if plot_over_time:
    print("Plot over time.")

  last_length = 0
  
  def init():
    global geometry_component, line_2D, lines_2D, lines_3D, line_comp, cbar, top_text, ax1, ax2, cmap, show_geometry, show_components, show_specified_fields, \
      specified_fields, solution_components, solution_name, solution_component, scaling_factors, min_x, max_x, min_y, max_y, min_z, max_z, compute_stress
      
    field_variable_names = py_reader.get_field_variable_names(data[0])
    
    # determine solution variable and component
    # if field names have been specified, only plot those
    if show_specified_fields:
    
      # get the specified fields as (field_variable_name, component_name) tuples
      specified_fields = []
      for specified_field_name in specified_field_names:
        for field_variable_name in field_variable_names:
          component_names = py_reader.get_component_names(data[0], field_variable_name)
          if field_variable_name == specified_field_name:
            if len(component_names) == 1:
              component_name = component_names[0]
              specified_fields.append((field_variable_name,component_name))
            else:
              for component_name in component_names:
                specified_fields.append((field_variable_name,component_name))
          elif specified_field_name in component_names:
            specified_fields.append((field_variable_name,specified_field_name))
      
      if len(specified_fields) == 0:
        print("\033[0;31mError! Given names {} do not match any existing field variable names or component names.\033[0m\n".format(specified_field_names))
        quit()
      
      n_components = 0
      
      # determine min/max values
      min_s = None
      max_s = None
      for (field_variable_name, component_name) in specified_fields:
        min_s_, max_s_ = py_reader.get_min_max(data, field_variable_name, component_name)
        print("\"{}.{}\" value range: [{}, {}]".format(field_variable_name,component_name, min_s_, max_s_))
      
        min_s = (min_s_ if min_s is None else min(min_s_, min_s))
        max_s = (max_s_ if max_s is None else max(max_s_, max_s))
         
    # if no field names have been specified, search for the solution variable and component
    else:
      solution_name = "solution"
      if "solution" not in field_variable_names:
        for field_variable_name in field_variable_names:
          if field_variable_name != "geometry":
            component_names = py_reader.get_component_names(data[0], field_variable_name)
            if len(component_names) == 1:
              solution_name = field_variable_name
              break
              
      component_names = py_reader.get_component_names(data[0], solution_name)
      solution_component = component_names[0]
      
      min_s, max_s = py_reader.get_min_max(data, solution_name, solution_component)
      print( "\"{}\" value range: [{}, {}]".format(solution_name, min_s, max_s))
        
      solution_components = py_reader.get_component_names(data[0], solution_name)
      n_components = len(solution_components)
      
    if min_s is None or np.isinf(min_s) or np.isnan(min_s):
      min_s = 0
    if max_s is None or np.isinf(max_s) or np.isnan(max_s):
      max_s = 0
    
    span_x = max_x - min_x
    span_y = max_y - min_y
    span_z = max_z - min_z
    
    # if the geometry is a line, do not show geometry plot
    if span_y == 0 and span_z == 0:
      show_geometry = False
    if (not show_geometry) and n_components > 1:
      if solution_components[1] != "1":
        if not show_specified_fields:
          show_components = True
    
    if span_x >= span_y and span_x >= span_z:
      geometry_component = "x"
    elif span_y >= span_x and span_y >= span_z:
      geometry_component = "y"
    else:
      geometry_component = "z"
    
    if plot_over_time:
      # x-axis is time
      min_x = data[0]['currentTime']
      max_x = data[-1]['currentTime']
    else:  
      # x-axis is geometry
      min_x, max_x = py_reader.get_min_max(data, "geometry", geometry_component)
    
    if show_geometry:
      print( "geometry bounding box: x:[{},{}], y:[{},{}], z:[{},{}]".format(min_x,max_x,min_y,max_y,min_z,max_z))
    
    # prepare plot
    if show_geometry:
      gs = gridspec.GridSpec(2,1,height_ratios=[4,1])
      ax2 = plt.subplot(gs[0], projection='3d')
      ax1 = plt.subplot(gs[1])
    elif show_components:
      gs = gridspec.GridSpec(2,1,height_ratios=[3,4])
      ax1 = plt.subplot(gs[0])  # main component
      ax3 = plt.subplot(gs[1])  # all other components
    else:
      ax1 = plt.gca()
    
    # prepare plot of specified fields
    if show_specified_fields:
      lines_2D = []
      for (field_variable_name, component_name) in specified_fields:
        if component_name == "0":
          label = field_variable_name
        else:
          label = component_name
        if len(data) > 100:
          line, = ax1.plot([], [], '+-', lw=2, label=label)
        else:
          line, = ax1.plot([], [], '-', lw=2, label=label)
        lines_2D.append(line)
      ax1.set_xlim(min_x, max_x)
      ax1.set_xlabel('t')
      plt.subplots_adjust(right=0.66, top=0.94, bottom=0.18)
      ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., frameon=False)
      
    else:
      # prepare main plot
      if data[0]["basisFunction"] == "Hermite" and data[0]["onlyNodalValues"] == False:  # for Hermite do not show nodes
        line_2D, = ax1.plot([], [], '-', color="b", lw=2, label=solution_components[0])
      elif data[0]["basisFunction"] == "Lagrange"  and data[0]["basisOrder"] == 2:   # for quadratic Lagrange also do not show nodes
        line_2D, = ax1.plot([], [], '-', color="b", lw=2, label=solution_components[0])
      elif len(data) > 100:
        line_2D, = ax1.plot([], [], '-', color="b", lw=2, label=solution_components[0])
      else:
        line_2D, = ax1.plot([], [], '+-', color="b", lw=2, label=solution_components[0])
      
      last_length = 0
    
      xlabel = geometry_component
      if plot_over_time:
        ax1.set_xlabel('t')
      else:
        ax1.set_xlabel(xlabel.upper())
      ax1.set_ylabel(solution_name)
      if solution_components[0] != "0":
        plt.subplots_adjust(right=0.66, top=0.94, bottom=0.18)
        ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., frameon=False) 
        plt.tight_layout()
        
    ax1.set_xlim(min_x, max_x)
    margin = abs(max_s - min_s) * 0.1
    ax1.set_ylim(min_s - margin, max_s + margin)
    top_text = ax1.text(0.5, 0.94, "", horizontalalignment='center', transform=ax1.transAxes, family='monospace')
    plt.grid(which='major')
      
    # prepare geometry plot
    if show_geometry:
      ax2.set_title("geometry")
      ax2.set_xlim(min_x, max_x)
      ax2.set_ylim(min_y, max_y)
      ax2.set_zlim(min_z, max_z)
      
      lines_3D = []
      # plot line segments with corresponding color
      xdata = py_reader.get_values(data[0], "geometry", "x")
      for i in range(len(xdata)):
        p, = ax2.plot(xdata[i:i+2], xdata[i:i+2], xdata[i:i+2], lw=3)
        lines_3D.append(p)
      ax2.set_xlabel('X')
      ax2.set_ylabel('Y')
      ax2.set_zlabel('Z')
        
      # manually create colorbar  [https://stackoverflow.com/questions/8342549/matplotlib-add-colorbar-to-a-sequence-of-line-plots]
      plt.sca(ax2)
      
      norm = mpl.colors.Normalize(vmin=min_s, vmax=max_s)
      cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.jet)
      cmap.set_array([])
      cbar = fig.colorbar(cmap)
    
    # prepare other components plot
    if show_components:
      # determine scaling factors and axis limits for plot 2
      scaling_factors = []
      line_comp = []
      min_value = 0
      max_value = 1
      for (j,component_name) in enumerate(solution_components):
        min_comp, max_comp = py_reader.get_min_max(data, solution_name, component_name)
        values_comp = py_reader.get_values(data[0], solution_name, component_name)
        
        v = max(abs(max_comp), abs(min_comp))
        if abs(v) < 1e-5:
          v = 1e-5
          
        scaling_factor = abs(1./v)
        scaling_factors.append(scaling_factor)
        
        #print "{}, scaling_factor: {}, min_comp: {}, max_comp: {}".format(j, scaling_factor, min_comp, max_comp)
      
        if j > 0:
          min_value = min(min_value, min_comp*scaling_factor)
          max_value = max(max_value, max_comp*scaling_factor)
        if len(data) > 100:
          line_plot, = ax3.plot([], [], '-',  lw=1, label=component_name)
        else:
          line_plot, = ax3.plot([], [], '+-', lw=1, label=component_name)
        line_comp.append(line_plot)
        
        #print "   min_value: {} -> {}, max_value: {} -> {}".format(min_comp*scaling_factor, min_value, max_comp*scaling_factor, max_value)
      
      # compute stress
      component_names = py_reader.get_component_names(data[0], solution_name)
      if "razumova/A_1" in component_names and "razumova/A_2" in component_names:
        compute_stress = True
        line_plot, = ax3.plot([], [], '-', color="red", lw=4, label="stress")
        line_comp.append(line_plot)
        min_A_2, max_A_2 = py_reader.get_min_max(data, solution_name, "razumova/A_2")
        min_stress = ((( (min_A_2/140)*0.05) - 0.000107)/0.0021)*0.840625
        max_stress = ((( (max_A_2/140)*0.05) - 0.000107)/0.0021)*0.840625
        print("min stress: {}, max stress: {}".format(min_stress,max_stress))
        
        v = max(abs(max_stress), abs(min_stress))
        if abs(v) < 1e-5:
          v = 1e-5
        scaling_factor = abs(1./v)
        scaling_factors.append(scaling_factor)
        
      ax3.set_xlim(min_x, max_x)
      if plot_over_time:
        ax3.set_xlabel('t')
      margin = abs(max_value - min_value) * 0.1
      ax3.set_ylim(min_value - margin, max_value + margin)
      ax3.set_ylabel('Other components')
      if len(solution_components) > 5:
        ncol = max(1,len(solution_components)/10)
        plt.subplots_adjust(right=0.66, top=0.94, bottom=0.18)
        ax3.legend(prop={'size': 6, }, ncol=ncol, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., frameon=False)
      else:
        plt.subplots_adjust(right=0.66, top=0.94, bottom=0.18)
        ax3.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., frameon=False)
    return top_text,

  def animate(i):
    global top_text, solution_name, solution_component, last_length, compute_stress
    
    ##################
    # 2D plot of main solution component
    # display data
    
    if show_specified_fields:
        for i,(field_variable_name, component_name) in enumerate(specified_fields):
          
          xdata = []
          sdata = []
          for d in data:
            solution_values = py_reader.get_values(d, field_variable_name, component_name)
            if solution_values is None:
              print("field_variable_name: {}, component_name: {}".format(field_variable_name, component_name))
              print("d: {}".format(d))
            xdata.append(d['currentTime'])
            sdata.append(solution_values[0])
          
          # refresh the line object that is the graph of the curve
          lines_2D[i].set_data(xdata,sdata)
      
    else:
      
      # plot over time instead of geometry (for cellml single instance)
      if plot_over_time:
        xdata = []
        sdata = []
        
        for d in data:
          solution_values = py_reader.get_values(d, solution_name, solution_component)
          
          xdata.append(d['currentTime'])
          sdata.append(solution_values[0])
      
      else:
        # plot over geometry
        xdata = py_reader.get_values(data[i], "geometry", geometry_component)
        sdata = py_reader.get_values(data[i], solution_name, solution_component)

        # handle Hermite that have derivative values saved
        if data[i]["basisFunction"] == "Hermite" and data[i]["onlyNodalValues"] == False:
          
          def hermite0(xi):
            return 1 - 3*xi*xi + 2*xi*xi*xi
            
          def hermite1(xi):
            return xi * (xi-1) * (xi-1)

          def hermite2(xi):
            return xi*xi * (3 - 2*xi)

          def hermite3(xi):
            return xi*xi * (xi-1)
          
          n_elements = data[i]["nElements"][0]
          
          n = 20
          new_xdata = np.zeros(n_elements*n)
          new_sdata = np.zeros(n_elements*n)
          
          #print("n entries: {}, new_xdata:{}".format(n_elements*n, new_xdata))
          #print("xdata: {}".format(xdata))
          
          for el_no in range(n_elements):
            c0 = sdata[2*el_no+0]
            c1 = sdata[2*el_no+1]
            c2 = sdata[2*el_no+2]
            c3 = sdata[2*el_no+3]
            
            #print("parsed coefficients: {} {} {} {}".format(c0,c1,c2,c3))
            
            for j in range(n):
              xi = float(j)/n
              x = xdata[2*el_no+0]*hermite0(xi) + xdata[2*el_no+1]*hermite1(xi) + xdata[2*el_no+2]*hermite2(xi) + xdata[2*el_no+3]*hermite3(xi)
              
              #print("xi={}, x={}".format(xi,x))
              
              new_xdata[el_no*n+j] = x
              new_sdata[el_no*n+j] = c0*hermite0(xi) + c1*hermite1(xi) + c2*hermite2(xi) + c3*hermite3(xi)
            
              #print("xi={}, s={:.2e}={:.2e}*{:.2e}+ {:.2e}*{:.2e}+ {:.2e}*{:.2e}+ {:.2e}*{:.2e}".format(xi,new_sdata[el_no*n+j],c0,hermite0(xi),c1,hermite1(xi),c2,hermite2(xi),c3,hermite3(xi)))
                          
          xdata = new_xdata
          sdata = new_sdata
            
        elif data[i]["basisFunction"] == "Lagrange"  and data[i]["basisOrder"] == 2:
            
          def q0(xi):
            return (2*xi - 1) * (xi-1)
          
          def q1(xi):
            return 4*(xi - xi*xi)
            
          def q2(xi):
            return 2*xi*xi - xi
            
          n_elements = data[i]["nElements"][0]
          
          n = 20
          new_xdata = np.zeros(n_elements*n)
          new_sdata = np.zeros(n_elements*n)
          
          #print("n entries: {}, new_xdata:{}".format(n_elements*n, new_xdata))
          #print("xdata: {}".format(xdata))
          
          for el_no in range(n_elements):
            c0 = sdata[2*el_no+0]
            c1 = sdata[2*el_no+1]
            c2 = sdata[2*el_no+2]
            
            #print("parsed coefficients: {} {} {}".format(c0,c1,c2))
            
            for j in range(n):
              xi = float(j)/n
              x = xdata[2*el_no+0]*q0(xi) + xdata[2*el_no+1]*q1(xi) + xdata[2*el_no+2]*q2(xi)
              
              #print("xi={}, x={} {}*{}+ {}*{}+ {}*{}".format(xi,x,xdata[2*el_no+0],q0(xi),xdata[2*el_no+1],q1(xi),xdata[2*el_no+2],q2(xi)))
              
              new_xdata[el_no*n+j] = x
              new_sdata[el_no*n+j] = c0*q0(xi) + c1*q1(xi) + c2*q2(xi)
            
          xdata = new_xdata
          sdata = new_sdata
          
      # refresh the line object that is the graph of the curve
      line_2D.set_data(xdata,sdata)
      
    ##################
    # 3D plot of geometry
    
    if show_geometry:
      # retrieve all values
      xdata = py_reader.get_values(data[0], "geometry", "x")
      ydata = py_reader.get_values(data[0], "geometry", "y")
      zdata = py_reader.get_values(data[0], "geometry", "z")
      
      min_s = min(sdata)
      max_s = max(sdata)
      
      # plot line segments with corresponding color
      for j in range(len(xdata)):
        normalized_value = (float)(sdata[j] - min_s) / (max_s - min_s)
        lines_3D[j].set_data([xdata[j:j+2], ydata[j:j+2]])
        lines_3D[j].set_3d_properties(zdata[j:j+2])
        lines_3D[j].set_color(plt.cm.jet(normalized_value))
        
    ##################
    # 2D plot of other components
    if show_components:
      for (j,component_name) in enumerate(solution_components):
        # do not plot main component
        if j == 0:
          continue
            
        # plot over time instead of geometry
        if plot_over_time:
          xdata = []
          data_comp = []
          if plot_over_time:
            for d in data:
              solution_values = py_reader.get_values(d, solution_name, component_name)
              
              xdata.append(d['currentTime'])
              data_comp.append(solution_values[0])
        else:
          data_comp = py_reader.get_values(data[i], solution_name, component_name)
          
        # refresh the line object that is the graph of the curve
        line_comp[j].set_data(xdata,np.array(data_comp)*scaling_factors[j])
      
      # compute stress
      if compute_stress:
        xdata = []
        data_comp = []
        for d in data:
          A_1 = py_reader.get_values(d, solution_name, "razumova/A_1")
          A_2 = py_reader.get_values(d, solution_name, "razumova/A_2")
          stress = ((( (A_1/140)*0.0+ (A_2/140)*0.05) - 0.000107)/0.0021)*0.840625
          
          xdata.append(d['currentTime'])
          data_comp.append(stress[0])
        
        # refresh the line object that is the graph of the curve
        line_comp[-1].set_data(xdata,np.array(data_comp)*scaling_factors[-1])
    
        
    # display timestep
    if not plot_over_time:
      if 'timeStepNo' in data[i]:
        timestep = data[i]['timeStepNo']
      if 'currentTime' in data[i]:
        current_time = data[i]['currentTime']
        
      max_timestep = len(data)-1
        
      t = "timestep {}/{}, t = {}".format(timestep, max_timestep, current_time)
      if last_length > len(t):
        t += " "*(last_length-len(t))
      last_length = len(t)
      top_text.set_text(t)
      
    return top_text,
    
  interval = 5000.0 / len(data)
        
  if len(data) == 1 or plot_over_time:
    init()
    animate(0)
    plt.savefig("fig.pdf")
    
  else:
    # create animation
    anim = animation.FuncAnimation(fig, animate, init_func=init,
               frames=len(data), interval=interval, blit=False)

    try:
      anim.save("anim.mp4")
    except:
      print("An error occured during the animation.")
    
    # create plot with first and last dataset
    # plot first dataset
    plt.clf()
    init()
    line_2D, = ax1.plot([], [], '+-', color=(1.0,0.9,0.8), lw=2, label="t=0")
    animate(0)
    
    # plot last dataset
    i = len(data)-1
    if 'currentTime' in data[i]:
      current_time = data[i]['currentTime']
    line_2D, = ax1.plot([], [], '+-', color=(1.0,0.9,0.8), lw=2, label="t={}".format(current_time))
    
    animate(i)
    
    max_timestep = len(data)-1
    if 'timeStepNo' in data[i]:
      timestep = data[i]['timeStepNo']
    top_text.set_text("timesteps 0 and {}".format(timestep))
    
    if show_geometry:
      ax2.set_title("geometry for t={}".format(current_time))
    
    plt.sca(ax1)
    #ax.add_line(line0)
    #ax.add_line(line1)
    plt.legend()
      
    plt.savefig("fig.pdf")
    
  if show_plot:
    plt.show()
  
####################
# 2D
if dimension == 2:
  
  field_variable_names = py_reader.get_field_variable_names(data[0])

  solution_name = "solution"
  if "solution" not in field_variable_names:
    for field_variable_name in field_variable_names:
      if field_variable_name != "geometry":
        component_names = py_reader.get_component_names(data[0], field_variable_name)
        if len(component_names) == 1:
          solution_name = field_variable_name
          break
          
  component_names = py_reader.get_component_names(data[0], solution_name)
  solution_component = component_names[0]
  
  # classical 2D scalar field variables, like in Laplace eq.
  if "displacements" not in field_variable_names:
    
    debug = False
    
    min_value, max_value = py_reader.get_min_max(data, solution_name, solution_component)
    min_x, max_x = py_reader.get_min_max(data, "geometry", "x")
    min_y, max_y = py_reader.get_min_max(data, "geometry", "y")
    
    print( "value range: [{}, {}]".format(min_value, max_value))
    
    # prepare plot
    fig = plt.figure(1)

    margin = abs(max_value - min_value) * 0.1
    ax = fig.add_subplot(111, projection='3d', xlim=(min_x, max_x), ylim=(min_y, max_y), zlim=(min_value-margin, max_value+margin))
    
    # surface = ax.plot_surface([], [], [], cmap=cm.coolwarm, linewidth=1,rstride=1,cstride=1) # needed? error with python3
    text = plt.figtext(0.15,0.85,"timestep",size=20)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    
    # create mesh
  
    if data[0]["basisFunction"] == "Lagrange":
      n_average_nodes_1D_per_element = data[0]["basisOrder"]
    
    elif data[0]["basisFunction"] == "Hermite":
      n_average_nodes_1D_per_element = 1
    
    if data[0]["meshType"] == "StructuredRegularFixed" or data[0]["meshType"] == "RegularFixed" or data[0]["meshType"] == "StructuredDeformable":
      
      if debug:
        print( "basisfunction: [{}], basisOrder: [{}]".format(data[0]["basisFunction"], data[0]["basisOrder"]))
      
      nEntries = []
      for i in range(dimension):
        if "nElements" in data[0]:
          nElements = data[0]["nElements"][i]
        else:
          nElements = data[0]["nElementsLocal"][i]
        nEntries.append(n_average_nodes_1D_per_element * nElements + 1)

      nEntries = nEntries[::-1]   # reverse list
      
      x_positions = py_reader.get_values(data[0], "geometry", "x")
      y_positions = py_reader.get_values(data[0], "geometry", "y")
      
      X = np.reshape(x_positions, nEntries)
      Y = np.reshape(y_positions, nEntries)
      
      if debug:
        print( "nEntries: ", nEntries)
        print( "x_positions: ", x_positions)
        print( "X: ",X)
        print( "y_positions: ", y_positions)
        print( "Y: ",Y)
      
      #print "x_positions shape: {}".format(len(x_positions))
      
    elif data[0]["meshType"] == "UnstructuredDeformable":
      if not data[0]["onlyNodalValues"]:
        print( "Error: onlyNodalValues is False, set to True in OutputWriter config!")
      
      x_positions = py_reader.get_values(data[0], "geometry", "x")
      y_positions = py_reader.get_values(data[0], "geometry", "y")
      X = x_positions
      Y = y_positions
      
      triangles = []
      for elemental_dofs in data[0]["elementalDofs"]:
        
        # for linear Lagrange and Hermite
        if n_average_nodes_1D_per_element == 1:
          triangles.append([elemental_dofs[0], elemental_dofs[1], elemental_dofs[3]])
          triangles.append([elemental_dofs[0], elemental_dofs[3], elemental_dofs[2]])
        else:  # for quadratic Lagrange
          triangles.append([elemental_dofs[0], elemental_dofs[1], elemental_dofs[4]])
          triangles.append([elemental_dofs[0], elemental_dofs[4], elemental_dofs[3]])
          triangles.append([elemental_dofs[4], elemental_dofs[1], elemental_dofs[2]])
          triangles.append([elemental_dofs[4], elemental_dofs[2], elemental_dofs[5]])
          triangles.append([elemental_dofs[4], elemental_dofs[5], elemental_dofs[6]])
          triangles.append([elemental_dofs[4], elemental_dofs[6], elemental_dofs[7]])
          triangles.append([elemental_dofs[4], elemental_dofs[6], elemental_dofs[3]])
          triangles.append([elemental_dofs[4], elemental_dofs[7], elemental_dofs[6]])
    
    def animate(i):
      ax.clear()
      
      # display data
      solution_shaped = py_reader.get_values(data[i], solution_name, solution_component)
      
      try:
        Z = np.reshape(solution_shaped, nEntries)
      except:
        Z = solution_shaped
      
      if debug:
        try:
          print( "x shape: {}, y shape: {}, z shape: {}".format(X.shape, Y.shape, Z.shape))
        except:
          pass
      
      # for unstructured grid use plot_trisurf
      if data[0]["meshType"] == "UnstructuredDeformable":
        plot = ax.plot_trisurf(X, Y, triangles, Z, cmap=cm.coolwarm, linewidth=1)
        
      # for structured grids use plot_surface
      else:
        plot = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=1,rstride=1,cstride=1)
        
      ax.set_zlim(min_value-margin, max_value+margin)
      ax.set_xlabel('X')
      ax.set_ylabel('Y')
      ax.set_zlabel('Z')
      
      # display timestep
      if 'timeStepNo' in data[i]:
        timestep = data[i]['timeStepNo']
        max_timestep = data[-1]['timeStepNo']
      if 'currentTime' in data[i]:
        current_time = data[i]['currentTime']
        
      if timestep == -1 or timestep == 0 or timestep == 1:
        text.set_text("t = {}".format(current_time))
      else:
        text.set_text("timestep {}/{}, t = {}".format(timestep, max_timestep, current_time))
      
      return plot,
      
    interval = 5000.0 / len(data)
          
    if len(data) == 1:
      animate(0)
      plt.savefig("fig.pdf")
      
    else:     # create animation
      anim = animation.FuncAnimation(fig, animate,
                 frames=len(data), interval=interval, blit=False)

      anim.save("anim.mp4")
      
      # create plot with first and last dataset
      fig2 = plt.figure(2,figsize=(5,10))
      ax2 = fig2.add_subplot(211, projection='3d', xlim=(min_x, max_x), ylim=(min_y, max_y), zlim=(min_value-margin, max_value+margin))
      
      ax1 = ax
      
      # plot first dataset
      ax = ax2
      plot0, = animate(0)
      
      # plot last dataset
      ax2 = fig2.add_subplot(212, projection='3d', xlim=(min_x, max_x), ylim=(min_y, max_y), zlim=(min_value-margin, max_value+margin))
    
      fig2.suptitle('first and last timesteps')
    
      i = len(data)-1
      if 'timeStepNo' in data[i]:
        timestep = data[i]['timeStepNo']
      if 'currentTime' in data[i]:
        current_time = data[i]['currentTime']
        
      ax = ax2
      plot1, = animate(i)
      
      ax = ax1
      
      #ax = plt.gca()
      #ax.add_line(line0)
      #ax.add_line(line1)
      plt.legend()
        
      plt.tight_layout()
      plt.savefig("fig.pdf")
      
    if show_plot:
      plt.show()
  
  elif "displacements" in field_variable_names:
  # continuum mechanics
    
    debug = False
    if debug:
      # print data
      import pprint 
      pp = pprint.PrettyPrinter()
      pp.pprint(data[0])
    
      print( "field_variable_names: ", field_variable_names)
    
    min_x_current, max_x_current = py_reader.get_min_max(data, "geometry", "x")
    min_y_current, max_y_current = py_reader.get_min_max(data, "geometry", "y")
    
    min_x_reference, max_x_reference = py_reader.get_min_max(data, "geometryReference", "0")
    min_y_reference, max_y_reference = py_reader.get_min_max(data, "geometryReference", "1")
    
    min_x = min(min_x_current, min_x_reference)
    min_y = min(min_y_current, min_y_reference)
    max_x = max(max_x_current, max_x_reference)
    max_y = max(max_y_current, max_y_reference)
  
    # prepare plot
    fig = plt.figure()

    if debug:
      print( "min_x: ", min_x)
      print( "min_y: ", min_y)
      print( "max_x: ", max_x)
      print( "max_y: ", max_y)

    # create plot with 30% margins
    margin_x = abs(max_x - min_x) * 0.3
    margin_y = abs(max_y - min_y) * 0.3
    ax = fig.add_subplot(111, xlim=(min_x-margin_x, max_x+margin_x), ylim=(min_y-margin_y, max_y+margin_y))
    text = plt.figtext(0.15,0.85,"timestep",size=20)
    
    ax.set_xlabel('x')
    ax.set_ylabel('y')
  
    def init():
      pass
        
    def animate(frame_no):
      ax.clear()
    
      dataset = data[frame_no]
      
      # create mesh / polygons
      polygons = []
      if dataset["meshType"] == "StructuredRegularFixed" or dataset["meshType"] == "StructuredDeformable":
        
        if debug:
          print( "basisfunction: [{}], basisOrder: [{}]".format(dataset["basisFunction"], dataset["basisOrder"]))
        
        if dataset["basisFunction"] == "Lagrange":
          nEntries = dimension * [0]
          for i in range(dimension):
            
            if "nElements" in dataset:
              nElements = dataset["nElements"][i]
            else:
              nElements = dataset["nElementsLocal"][i]
            nEntries[i] = dataset["basisOrder"] * nElements + 1
            
        elif dataset["basisFunction"] == "Hermite":
          nEntries = dimension * [0]
          for i in range(dimension):
            if "nElements" in dataset:
              nElements = dataset["nElements"][i]
            else:
              nElements = dataset["nElementsLocal"][i]
            nEntries[i] = nElements + 1
        
        nEntries = nEntries[::-1]   # reverse list
      
      def create_mesh(x_positions, y_positions, nEntries, **kwargs):
        """ create Polygon objects (quads) that can be plotted """
        polygons = []
        
        X = np.reshape(list(x_positions), nEntries)
        Y = np.reshape(list(y_positions), nEntries)
        
        if debug:
          print( "nEntries: ", nEntries)
          print( "x_positions: ", x_positions)
          print( "X: ",X)
          print( "y_positions: ", y_positions)
          print( "Y: ",Y)
          print( nEntries)
        
        # loop over elements
        for ely in range(nEntries[0]-1):
          for elx in range(nEntries[1]-1):
            point0 = np.array([X[ely][elx], Y[ely][elx]])
            point1 = np.array([X[ely][elx+1], Y[ely][elx+1]])
            point2 = np.array([X[ely+1][elx], Y[ely+1][elx]])
            point3 = np.array([X[ely+1][elx+1], Y[ely+1][elx+1]])
            
            if debug:
              print( "polygon (0,1,2,3) = ({},{},{},{})".format(point0, point1, point2, point3))
        
            polygon = Polygon([point0, point1, point3, point2], **kwargs)  # dummy data for xs,ys
            polygons.append(polygon)
            
        return polygons
      
      # parse values of reference and current configuration
      x_positions_current = py_reader.get_values(dataset, "geometry", "x")
      y_positions_current = py_reader.get_values(dataset, "geometry", "y")
      
      x_positions_reference = py_reader.get_values(dataset, "geometryReference", "0")
      y_positions_reference = py_reader.get_values(dataset, "geometryReference", "1")
        
      # create meshes for current and reference configuration
      polygons_current = []
      polygons_current = create_mesh(x_positions_current, y_positions_current, nEntries, fill=False, edgecolor=(0.2, 0.2, 0.8), linewidth=3, label="current configuration")
      polygons_reference = []
      polygons_reference = create_mesh(x_positions_reference, y_positions_reference, nEntries, fill=False, edgecolor=(0.8,0.8,0.8), linewidth=3, label="reference configuration")
      
      
      # add all polygons to axis
      for polygon in polygons_reference + polygons_current:
        ax.add_patch(polygon)
        
      ax.set_aspect('equal','datalim')

      # add legend, every legend entry (current configuration or reference configuration) only once
      handles, labels = ax.get_legend_handles_labels()
      new_handles = []
      new_labels = []
      for (handle,label) in zip(handles,labels):
        if label not in new_labels:
          new_handles.append(handle)
          new_labels.append(label)
      ax.legend(new_handles, new_labels, loc='best')
      
     
      # display timestep
      if 'timeStepNo' in dataset:
        timestep = dataset['timeStepNo']
        max_timestep = data[-1]['timeStepNo']
      if 'currentTime' in dataset:
        current_time = dataset['currentTime']
        
      t = "{}. timestep {}/{}, t = {}".format(frame_no, timestep, max_timestep, current_time)
      text.set_text(t)
      #text = plt.figtext(0.15, 0.85, t, size=20)
    
    interval = 5000.0 / len(data)
          
    if len(data) == 1:
      init()
      animate(0)
      plt.savefig("fig.pdf")
      
      if show_plot:
        plt.show()
      
    else:  # create animation
      
      anim = animation.FuncAnimation(fig, animate, init_func=init,
                 frames=len(data), interval=interval, blit=False)
        
      if show_plot:
        plt.show()
      else:
        anim.save("anim.mp4")
        plt.savefig("fig.pdf")
      
sys.exit(0)
