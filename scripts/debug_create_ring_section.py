#!/usr/bin/env ../../dependencies/python/install/bin/python3 
# -*- coding: utf-8 -*-
#
# This script calls the create_ring_section function from geometry_manipulation/stl_create_rings.py

import sys,os
import pickle
import stl_create_rings

import datetime
now = datetime.datetime.now()
print(" ======= debug_create_ring_section.py =======") 
print(now.strftime("%d/%m/%Y %H:%M:%S"))

# parse command line arguments
input_filename = "/store/software/opendihu/testing/system_testing/tests/fibers/meshes/biceps_full.stl"

# parameters
# for initial mesh
bottom_clip = -600.0     # bottom clip plane of triangles that will not be considered
top_clip = -290.0        # top clip plane

# for repaired mesh ("biceps_full.stl")
bottom_clip = 37.0
top_clip = 300.0

start_point = [70,  148,  48]
end_point = [66,  143,   51]
z_value = 50
n_points = 10

stl_create_rings.create_ring_section(input_filename, start_point, end_point, z_value, n_points)
