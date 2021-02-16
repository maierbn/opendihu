#!/usr/bin/env python
# -*- coding: utf-8 -*-

import datetime
now = datetime.datetime.now()
print(" ======= plot_mesh_quality.py =======") 
print(now.strftime("%d/%m/%Y %H:%M:%S"))

import sys, os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv

# load data
if len(sys.argv) == 2:
  filename = sys.argv[1]
else:
  print("usage: <filename>")
  #quit()

data = np.genfromtxt(filename, delimiter=";")
if False:
  data = np.array([
    [0,0,961,0.02190327871339396,0.1703234817249599], # 31x31fibers_l1_m1_4In_neumann.no_boundary.bin
    [0,0,961,0.021851135102074922,0.17050822677704625], # 31x31fibers_l1_m1_3In_neumann.no_boundary.bin
    [0,0,961,0.021813104172269338,0.17081704640689485], # 31x31fibers_l1_m1_2In_neumann.no_boundary.bin
    [0,0,961,0.022020439851158917,0.17388858740852675], # 31x31fibers_l1_m1_3In_dirichlet.no_boundary.bin
    [0,0,961,0.021958865533939623,0.1749230361859109], # 31x31fibers_l1_m1_2In_dirichlet.no_boundary.bin
    [0,0,961,0.025242259874085154,0.17506040702962086], # 31x31fibers_l1_m1_1IG_dirichlet.no_boundary.bin
    [0,0,961,0.022535321251316052,0.17586747159571584], # 31x31fibers_l1_m1_2IG_dirichlet.no_boundary.bin
    [0,0,961,0.02166185838487164,0.1759734570956548],   # 31x31fibers_l1_m1_1In_dirichlet.no_boundary.bin
  ])
  #"N4s", "N3s", "N2s", "D3s", "D2s", "D1g", "D2g", "D1s"
  
  data = np.array([
    [0,0,225,0.008933518485395723,0.10890216755506996],  # 15x15fibers_l1_m0_q2In_dirichlet.no_boundary.bin
    [0,0,225,0.008851797237756283,0.10935338702608419],  # 15x15fibers_l1_m0_q1In_dirichlet.no_boundary.bin
    [0,0,225,0.009155474976393472,0.10946341348996944],  # 15x15fibers_l1_m0_q3In_dirichlet.no_boundary.bin
    [0,0,225,0.008644140131063216,0.11012143439787514],  # 15x15fibers_l1_m0_1In_dirichlet.no_boundary.bin
    [0,0,225,0.008831287781090552,0.11025683526901013],  # 15x15fibers_l1_m0_3In_dirichlet.no_boundary.bin
    [0,0,225,0.008731135691359585,0.11034083734498665],  # 15x15fibers_l1_m0_2In_dirichlet.no_boundary.bin
    [0,0,225,0.06649357250892905,0.11567386103250858],  # 15x15fibers_l1_m0_q1In_neumann.no_boundary.bin
    [0,0,225,0.04088688335413133,0.17461451733968913],  # 15x15fibers_l1_m0_1nn_dirichlet.no_boundary.bin
    [0,0,225,0.04153353427598159,0.17607887140149636],  # 15x15fibers_l1_m0_3nn_dirichlet.no_boundary.bin
    [0,0,225,0.041404557510790994,0.17616284009062344],  # 15x15fibers_l1_m0_q2nn_dirichlet.no_boundary.bin
    [0,0,225,0.04106884743151952,0.17672550916210167],  # 15x15fibers_l1_m0_2nn_dirichlet.no_boundary.bin
    [0,0,225,1.7090131999185547,0.1785289173932961],  # 15x15fibers_l1_m0_1In_neumann.no_boundary.bin
    [0,0,225,0.04499659041048243,0.1797487435435816],  # 15x15fibers_l1_m0_q1nn_dirichlet.no_boundary.bin
    [0,0,225,0.4850027161025263,0.18075899915948593],  # 15x15fibers_l1_m0_3In_neumann.no_boundary.bin
    [0,0,289,0.01199135064045919,0.1820827679560739],  # 17x17fibers_l1_m0_q2In_dirichlet.bin
    [0,0,289,0.011882573915533193,0.182155291891949],  # 17x17fibers_l1_m0_q1In_dirichlet.bin
    [0,0,289,0.011879313079233452,0.18280965305112426],  # 17x17fibers_l1_m0_3In_dirichlet.bin
    [0,0,289,0.01163668497348027,0.1836279129147605],  # 17x17fibers_l1_m0_2In_dirichlet.bin
    [0,0,289,0.01214339740605334,0.18384885155415567],  # 17x17fibers_l1_m0_q3In_dirichlet.bin
  ])
# columns:
# 0 scenario name; filename; n valid fibers; variance element lengths; variance angles
# 1 filename; 
# 2 n valid fibers; 
# 3 variance element lengths; 
# 4 variance angles; 

labels = [
  "qD2s", "qD1s", "qD3s", "lD1s", "lD2s", "lD3s", "qN1s"
]

plt.rcParams.update({'font.size': 18})

tlist = np.arange(len(labels))  # the label locations
width = 0.7  # the width of the bars

#datalist = [data[5,4], data[4,4], data[3,4]]
datalist = data[0:len(labels),4]

textcolor = "#800000"
textcolor = "k"
color = "#ffcc00"
fig, ax = plt.subplots(figsize=(4.5,5))
rects1 = ax.barh(tlist, datalist, width, color=color)
for t,d in zip(tlist,datalist):
  ax.text(d-0.01,t-0.05,"{:.4f}".format(d),color=textcolor,verticalalignment='center',horizontalalignment='right')

ax.set_xlabel('Variance of angles [rad${}^2$]')
#ax.set_ylabel(r'')
ax.set_yticks(tlist)
ax.set_yticklabels(labels)
plt.grid(which='major', axis='x')
#ax.set_ylim(0,0.8)

#plt.legend([rects1[0]], ['Standard deviation of \nrelative element lengths', 'Duration [s]'], fontsize=16, loc='upper center')
# Add some text for labels, title and custom x-axis tick labels, etc.

fig.tight_layout()
plt.savefig("mesh_quality_options.pdf")
plt.show()
