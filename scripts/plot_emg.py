#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script to visualize EMG output files
# arguments: <filename>
#

import sys, os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv
import os.path

search_paths = [".", "build_release/out/biceps", "build_release/out", "build_release"]

emg_filename = "electrodes.csv"
stimulation_times_filename = "stimulation.log"

# parse file
if len(sys.argv) > 2:
  emg_filename = sys.argv[1]
  stimulation_times_filename = sys.argv[2]

elif len(sys.argv) > 1:
  emg_filename = sys.argv[1]
  
else:
  print("usage: {} <electrodes filename> <stimulation_times_filename> [<tstart>, <tend>]\n"
    "For example:\n"
    "  plot_emg.py electrodes.csv stimulation.log\n"
    "If no filename is given, search in subdirectories {}\n\n".format(sys.argv[0], search_paths))
  
  # get all input data in current directory
  for path in search_paths:
    try:
      files_at_path = os.listdir(path)

      if "electrodes.csv" in files_at_path:
        emg_filename = os.path.join(path, "electrodes.csv")
      if "stimulation.log" in files_at_path:
        stimulation_times_filename = os.path.join(path, "stimulation.log")
    except:
      pass

tstart = 0
tend = None
if len(sys.argv) == 5:
  tstart = (float)(sys.argv[3])
  tend = (float)(sys.argv[4])

# get filesize
try:
  filesize = os.path.getsize(emg_filename)
except:
  print("Could not open emg file \"{}\".".format(emg_filename))
  quit()

print("Input files:\n  EMG:         {}\n  stimulation: {}".format(emg_filename, stimulation_times_filename))
if tend is None:
  print("Plot whole time range. Use for example `plot_emg.py {} {} 1 100` to only plot the range 1-100 ms.".format(emg_filename, stimulation_times_filename))
else:
  print("Plot time range [{}, {}] ms.".format(tstart, tend))

# CC BY-SA 3.0 by Sridhar Ratnakumar https://stackoverflow.com/a/1094933/10290071
def sizeof_fmt(num, suffix='B'):
  for unit in ['','Ki','Mi','Gi','Ti','Pi','Ei','Zi']:
    if abs(num) < 1024.0:
      return "%3.1f %s%s" % (num, unit, suffix)
    num /= 1024.0
  return "%.1f%s%s" % (num, 'Yi', suffix)
print("Input file \"{}\" is {}".format(emg_filename, sizeof_fmt(filesize)), end="")

# depending on the settings "enableGeometryInCsvFile", there are two possible file formats of the csv file:
#"enableGeometryInCsvFile": True
# file starts with:
#|  #timestamp;t;n_points;p0_x;p0_y;p0_z;p1_x;p1_y;p1_z;p2_x;  (...)   p171_z;p0_value;p1_value;  (...)
#|  2020/9/28 09:16:40;0;172;9.01393;15.7361;-54.7243;         (...)

#"enableGeometryInCsvFile": True
# file starts with:
#|  #electrode positions (x0,y0,z0,x1,y1,z1,...);
#|  #; ;9.92639;13.9759;-54.8722;                   (...)
#|  #timestamp;t;n_points;p0_value;p1_value;        (...)
#|  2020/9/29 10:08:48;0;384;0.0030616;0.00300943;  (...)


# parse position data
position_data = None
file_contains_transient_position_data = True
electrode_data_begin = 3

# open csv file and check if the position data is contained only as comment before the actual data
with open(emg_filename, "r") as f:
  line = f.readline()
  if "#electrode positions" in line:
    line = f.readline()
    position_data = list(map(float,line.split(";")[2:]))
    file_contains_transient_position_data = False

# parse all csv values
data = np.genfromtxt(emg_filename, delimiter=';')
t_list = data[:,1]

# extract time range
index_start = 0
index_end = len(t_list)
for i,t in enumerate(t_list):
  if t < tstart:
    index_start = i
  if tend and t > tend:
    index_end = i
    break
  
t_list = t_list[index_start:index_end]
data = data[index_start:index_end,:]

# determine the first column that contains EMG values
n_points = (int)(data[0][2])
if file_contains_transient_position_data:
  electrode_data_begin = 3 + n_points*3
else:
  electrode_data_begin = 3

# determine number of electrodes in x direction
# extract all z positions (along muscle)
if file_contains_transient_position_data:
  position_data = data[0][3:3+3*n_points]
electrode_positions = [np.array(position_data[3*i:3*i+3]) for i in range(n_points)]
z_positions = np.array([position_data[3*i+2] for i in range(n_points)])

# compute differences in z position between neighbour points, this will have a jump (high increment) where a new line starts
differences = z_positions[1:] - z_positions[0:-1]

# determine the indices where the increment is higher than 3 times the mean 
jumps = np.where(differences > np.mean(differences)*3)[0]

# number of points in x/y-direction (across muscle) is position of first jump + 1
n_points_xy = jumps[0] + 1
n_points_z = (int)(n_points / n_points_xy)

# plot electrode positions for debugging
if False:
  from mpl_toolkits.mplot3d import Axes3D
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')

  xpos = [electrode_positions[i][0] for i in range(n_points)]
  ypos = [electrode_positions[i][1] for i in range(n_points)]
  zpos = [electrode_positions[i][2] for i in range(n_points)]
  colors = [(i/n_points,i/n_points,i/n_points) for i in range(n_points)]
  ax.scatter(xpos,ypos,zpos,c=colors)
  plt.show()

# compute average interelectrode distances
# across muscle (xy)
distances_xy = []
for j in range(n_points_z):
  for i in range(n_points_xy-1):
    point_index = j*n_points_xy + i
    
    p0 = electrode_positions[point_index]
    p1 = electrode_positions[point_index+1]
    
    distance = np.linalg.norm(p0-p1)
    distances_xy.append(distance)

# along muscle (z)
distances_z = []
for j in range(n_points_z-1):
  for i in range(n_points_xy):
    point_index = j*n_points_xy + i
    
    p0 = electrode_positions[point_index]
    p1 = electrode_positions[point_index+n_points_xy]
    
    distance = np.linalg.norm(p0-p1)
    distances_z.append(distance)

# compute inter eletrcode distances in cm
inter_electrode_distance_xy = np.mean(distances_xy)
inter_electrode_distance_z = np.mean(distances_z)

print(" and contains {} x {} = {} points".format(n_points_xy, n_points_z, n_points))
print("Inter-electrode distance: {:.2f} mm x {:.2f} mm".format(10*inter_electrode_distance_xy, 10*inter_electrode_distance_z))

# create grid plot
n_plots_x = n_points_xy
n_plots_y = n_points_z
fig = plt.figure(figsize=(5,10))

# set global parameters for font sizes
plt.rcParams.update({'font.size': 14})
#plt.rcParams['lines.linewidth'] = 1
#plt.rcParams['lines.markersize'] = 8

# labels for entire plot
ax = fig.add_subplot(111)
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
plt.grid(False)
plt.xlabel("cross-fiber direction")
plt.ylabel("fiber direction")  
plt.title("Surface EMG for {} x {} electrodes\n time range: [{}, {}] ms".format(n_points_xy, n_points_z, t_list[0], t_list[-1]))

# determine minimum and maximum overall values
emg_data = data[:,electrode_data_begin:]
minimum_value = np.min(emg_data)
maximum_value = np.max(emg_data)
print("EMG value range [mV]: [{},{}]".format(minimum_value, maximum_value))
print("time range [ms]: [{},{}]".format(t_list[0], t_list[-1]))
sampling_frequency = 1000*len(t_list) / (t_list[-1]-t_list[0])
print("sampling frequency [Hz]: {:.0f}".format(sampling_frequency))

# loop over sub plots
for j in range(n_plots_y):
  for i in range(n_plots_x):
    index = j*n_plots_x + i
    
    emg_list = emg_data[:,index]
    
    ax = fig.add_subplot(n_plots_y,n_plots_x,index+1)   # nrows, ncols, index, plots start top left and increase to the right
    ax.plot(t_list, emg_list)
    ax.set_ylim(minimum_value,maximum_value)
    #ax.axis('off')
    ax.set_xticks([])
    ax.set_yticks([])

plt.axis('off')
plt.subplots_adjust(hspace=0, wspace=0)
plt.subplots_adjust(left=0.17,right=0.928, bottom = 0.077, top = 0.929)
plt.savefig("emg_plot.pdf")
print("Created file \"emg_plot.pdf\".")

plt.ion()
plt.show()
plt.draw()
plt.pause(5)
plt.ioff()

# -------------------------
# load stimulation times
fiber_times = {}
fiber_mu_nos = {}

mu_times = {}
end_time = 0

# load data
try:
  with open(stimulation_times_filename) as f:
    for line in f:
      if '#' in line:
        continue
      values = line.split(";")
      mu_no = int(values[0])
      fiber_no = int(values[1])
      times = [float(v)/1000.0 for v in values[2:-1]]
      
      fiber_mu_nos[fiber_no] = mu_no
      fiber_times[fiber_no] = times
      
      if not mu_no in mu_times:
        mu_times[mu_no] = set()
      
      for item in times:
        if item > end_time:
          end_time = item
        mu_times[mu_no].add(item)

  print("Stimulation file \"{}\" contains end time {:.3f} s, {} motor units, {} fibers".format(stimulation_times_filename, end_time, len(mu_times), len(fiber_times)))

  max_mu = max(mu_times)
except:
  print("Stimulation file \"{}\" does not exit".format(stimulation_times_filename))
  max_mu = 1

# ------------------
# create animation
from matplotlib import animation
from matplotlib import gridspec

plt.rcParams.update({'font.size': 10})
fig = plt.figure(figsize=(8,10))
gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])

handles = []
time_line_handle = None
text_handle = None
colorbar_initialized = False

def init():
  global time_line_handle, text_handle, colorbar_initialized
  time_index = 0
  
  # top plot:
  # plot electrode positions as '+'
  # loop over electrodes
  for j in range(n_points_z):
    for i in range(n_points_xy):
      electrode_index = j*n_points_xy + i
  
      emg_data[time_index,electrode_index]
      handle = ax1.plot(i,j,'+',color=(0.5,0.5,0.5))
  
  # plot color
  data = np.transpose(np.reshape(emg_data[time_index,:],(n_points_xy,n_points_z)))
  handle = ax1.pcolormesh(data, cmap='RdBu', vmin=minimum_value, vmax=maximum_value)
  handles.append(handle)
  if not colorbar_initialized:
    fig.colorbar(handle, ax=ax1)
    colorbar_initialized = True
  ax1.set_aspect('equal')
  
  # left text
  plt.figtext(0.05, 0.95, "{} x {} = {} electrodes".format(n_points_xy, n_points_z, n_points), horizontalalignment='left')
  plt.figtext(0.05, 0.92, "IED: {:.1f} mm x {:.1f} mm".format(10*inter_electrode_distance_xy, 10*inter_electrode_distance_z), horizontalalignment='left')
  plt.figtext(0.05, 0.89, "t: {} ms".format(t_list[-1]-t_list[0]), horizontalalignment='left')
  plt.figtext(0.05, 0.86, "sampling frequency: {:.0f} Hz".format(sampling_frequency), horizontalalignment='left')
  plt.figtext(0.05, 0.83, "{} motor units".format(len(mu_times.items())), horizontalalignment='left')
  
  # bottom plot:
  # vertical line
  time_line_handle, = ax2.plot([time_index,time_index], [0,1], color='b')
  ax2.set_xlim(t_list[0]/1000, t_list[-1]/1000)
  ax2.set_ylim(0,max_mu)
  
  # text
  text_handle = ax2.text(0.05, 1.05, "", horizontalalignment='left', transform=ax2.transAxes, family='monospace')
  
  # firing times for fibers on mus
  if mu_times:

    # plot line for each MU
    for (mu_no,times) in mu_times.items():
      ax2.plot([0,end_time],[mu_no,mu_no],c=(0.8,0.8,0.8))

    for (fiber_no,times) in fiber_times.items():
      
      mu_no = fiber_mu_nos[fiber_no]
      # plot stimulation times as cross
      ax2.plot(list(times),[mu_no for _ in times],'+')
      
    ax2.set_ylabel('motor unit no [-]')
  ax2.set_xlabel('time [s]')
  
def animate(time_index):
  global time_line_handle, text_handle
  
  time_index = time_index*frame_stride
  
  if time_index % 100 == 0:
    print("{:.1f} %".format(100*time_index / len(t_list)))
  
  # top plot:
  data = emg_data[time_index,:]
  handles[0].set_array(data.flatten())
  
  # bottom plot:
  t = t_list[time_index]/1000
    
  time_line_handle.set_data([t,t], [0,max_mu])
  text_handle.set_text("{} s".format(t))
  

# compute timing values for the animation
slowdown_factor = 10   # how much slower the video will be than the actual simulation
if len(t_list) > 40000:  # more than 20s
  slowdown_factor = 2
  
target_fps = 20       # the framerate of the video

# duration of the video
duration = (t_list[-1] - t_list[0])/1000. * slowdown_factor
frame_stride = (int)(len(t_list) / (target_fps*duration))   # n_frames / duration = 20 => (len(t_list) // frame_stride) / ((t_list[-1] - t_list[0])/1000. * slowdown_factor) = 20 =>  frame_stride = len(t_list) // (20*((t_list[-1] - t_list[0])/1000. * slowdown_factor))
n_frames = len(t_list) // frame_stride
interval = (int)(duration / n_frames * 1000)

anim = animation.FuncAnimation(fig, animate, init_func=init, 
                               frames=n_frames, interval=interval, blit=False, repeat=False)

print("Saving animation with {} frames from {} data points, animation length {} s, data length {} s, {:.1f} frames/s".format(n_frames, len(t_list), duration, (t_list[-1] - t_list[0])/1000., n_frames/duration))

plt.tight_layout()

try:
  anim.save("anim.mp4")
  print("Saved \"anim.mp4\".")
except:
  print("An error occured during the animation.")

plt.show()

