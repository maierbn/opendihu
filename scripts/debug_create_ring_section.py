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

start_point=(66.55195617675781, 149.58843994140625, 60.0)
end_point=(76.11906433105469, 154.0915985107422, 60.0)
z_value=60.0
n_points=9

start_point=(76.11097717285156, 139.17101287841797, 61.89999771118164)
end_point=(85.38225555419922, 143.78616333007812, 60.000003814697266)
z_value=61.9
n_points=9


start_point=(53.1998,170.785,146)
end_point=(64.406,187.789,146)
z_value=72.0
n_points=9
result = stl_create_rings.create_ring_section(input_filename, start_point, end_point, z_value, n_points)

print("n points: {}, should be {}".format(len(result), n_points))
