# This is a helper script that sets a lot of the internal variables which are all defined in variables.py

import numpy as np
import scipy
import pickle
import sys,os
import struct
import argparse
import random
import time
sys.path.insert(0, '..')
import variables    # file variables.py

# parse arguments
rank_no = int(sys.argv[-2])
n_ranks = int(sys.argv[-1])

variables.n_subdomains = variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z

if variables.n_subdomains != n_ranks:
  print("\n\n\033[0;31mError! Number of ranks {} does not match given partitioning {} x {} x {} = {}.\033[0m\n\n".format(n_ranks, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z))
  quit()


def n_fibers_in_subdomain_x(_):
    return variables.n_fibers_x
def n_fibers_in_subdomain_y(_):
    return variables.n_fibers_y
def get_fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y):
    return fiber_in_subdomain_coordinate_x + variables.n_fibers_x * fiber_in_subdomain_coordinate_y

variables.n_subdomains_xy = variables.n_subdomains_x * variables.n_subdomains_y
variables.n_fibers_total = variables.n_fibers_x * variables.n_fibers_y

#############################
# create the partitioning using the script in create_partitioned_meshes_for_settings.py
# 3D meshes are created internally
# however, we have to create the fibers on our own

# fiber directions
fiber_meshes = {}
fiber_mesh_names = []

### muscle 1
for j in range(variables.n_fibers_y):
  for i in range(variables.n_fibers_x):
    fiber_no = j*variables.n_fibers_x + i
    
    # determine start position of fiber in (x,y)-plane
    x = 0 + i / (variables.n_fibers_x - 1) * variables.muscle_extent[0]
    y = 0 + j / (variables.n_fibers_y - 1) * variables.muscle_extent[1]

    # loop over points of a single fiber
    node_positions = []
    for k in range(variables.n_points_whole_fiber):
      x_pos = x
      y_pos = y
      z_pos = variables.muscle_offset[2] + k / (variables.n_points_whole_fiber - 1) * variables.muscle_extent[2]
      node_positions.append([x_pos,y_pos,z_pos])
    
    mesh_name = "muscle1_fiber{}".format(fiber_no)
    fiber_mesh_names.append(mesh_name)
    
    fiber_meshes[mesh_name] = {
      "nodePositions": node_positions,
      "nElements": [variables.n_points_whole_fiber - 1],
      "inputMeshIsGlobal": True,
      "nRanks": [n_ranks],
    }

#### set output writer
variables.output_writer_fibers_muscle1 = []
variables.output_writer_emg_muscle1 = []

### muscle 2
for j in range(variables.n_fibers_y):
  for i in range(variables.n_fibers_x):
    fiber_no = j*variables.n_fibers_x + i
    
    # determine start position of fiber in (x,y)-plane
    x = 0 + i / (variables.n_fibers_x - 1) * variables.muscle_extent[0]
    y = 0 + j / (variables.n_fibers_y - 1) * variables.muscle_extent[1]

    # loop over points of a single fiber
    node_positions = []
    for k in range(variables.n_points_whole_fiber):
      x_pos = x
      y_pos = y
      z_pos = variables.muscle_offset[2] + k / (variables.n_points_whole_fiber - 1) * variables.muscle_extent[2]
      node_positions.append([x_pos,y_pos,z_pos])
    
    mesh_name = "muscle2_fiber{}".format(fiber_no)
    fiber_mesh_names.append(mesh_name)
    
    fiber_meshes[mesh_name] = {
      "nodePositions": node_positions,
      "nElements": [variables.n_points_whole_fiber - 1],
      "inputMeshIsGlobal": True,
      "nRanks": [n_ranks],
    }

#### set output writer
variables.output_writer_fibers_muscle2 = []
variables.output_writer_emg_muscle2 = []

# set variable mappings for cellml model
if "hodgkin_huxley" in variables.cellml_file and "hodgkin_huxley-razumova" not in variables.cellml_file:
  # parameters: I_stim
  variables.muscle1_mappings = {
    ("parameter", 0):           ("constant", "membrane/i_Stim"),      # parameter 0 is constant 2 = I_stim
    ("connectorSlot", 0): ("state", "membrane/V"),              # expose state 0 = Vm to the operator splitting
  }
  variables.parameters_initial_values = [0.0]                         # initial value for stimulation current
  variables.nodal_stimulation_current = 40.                           # not used
  variables.vm_value_stimulated = 20.                                 # to which value of Vm the stimulated node should be set (option "valueForStimulatedPoint" of FastMonodomainSolver)

elif "shorten" in variables.cellml_file:
  # parameters: stimulation current I_stim, fiber stretch λ
  variables.muscle1_mappings = {
    ("parameter", 0):           ("algebraic", "wal_environment/I_HH"), # parameter is algebraic 32
    ("parameter", 1):           ("constant", "razumova/L_x"),             # parameter is constant 65, fiber stretch λ, this indicates how much the fiber has stretched, 1 means no extension
    ("connectorSlot", 0): ("state", "wal_environment/vS"),          # expose state 0 = Vm to the operator splitting
  }
  variables.parameters_initial_values = [0.0, 1.0]                        # stimulation current I_stim, fiber stretch λ
  variables.nodal_stimulation_current = 1200.                             # not used
  variables.vm_value_stimulated = 40.                                 # to which value of Vm the stimulated node should be set (option "valueForStimulatedPoint" of FastMonodomainSolver)
  
elif "slow_TK_2014" in variables.cellml_file:   # this is (3a, "MultiPhysStrain", old tomo mechanics) in OpenCMISS
  # parameters: I_stim, fiber stretch λ
  variables.muscle1_mappings = {
    ("parameter", 0):           ("constant", "wal_environment/I_HH"), # parameter 0 is constant 54 = I_stim
    ("parameter", 1):           ("constant", "razumova/L_S"),         # parameter 1 is constant 67 = fiber stretch λ
    ("connectorSlot","vm"):     "wal_environment/vS",                 # expose state 0 = Vm to the operator splitting
    ("connectorSlot", "stress"):"razumova/stress",                    # expose algebraic 12 = γ to the operator splitting
  }
  variables.parameters_initial_values = [0.0, 1.0]                    # wal_environment/I_HH = I_stim, razumova/L_S = λ
  variables.nodal_stimulation_current = 40.                           # not used
  variables.vm_value_stimulated = 40.                                 # to which value of Vm the stimulated node should be set (option "valueForStimulatedPoint" of FastMonodomainSolver)
  
elif "Aliev_Panfilov_Razumova_2016_08_22" in variables.cellml_file :   # this is (3, "MultiPhysStrain", numerically more stable) in OpenCMISS, this only computes A1,A2,x1,x2 not the stress
  # parameters: I_stim, fiber stretch λ, fiber contraction velocity \dot{λ}
  variables.muscle1_mappings = {
    ("parameter", 0):           ("constant", "Aliev_Panfilov/I_HH"),  # parameter 0 is constant 0 = I_stim
    ("parameter", 1):           ("constant", "Razumova/l_hs"),        # parameter 1 is constant 8 = fiber stretch λ
    ("parameter", 2):           ("constant", "Razumova/velo"),        # parameter 2 is constant 9 = fiber contraction velocity \dot{λ}
    ("connectorSlot", 0): ("state", "Aliev_Panfilov/V_m"),      # expose state 0 = Vm to the operator splitting
    ("connectorSlot", 1): ("algebraic", "Razumova/sigma"),   # expose algebraic 0 = γ to the operator splitting
  }
  variables.parameters_initial_values = [0, 1, 0]                     # Aliev_Panfilov/I_HH = I_stim, Razumova/l_hs = λ, Razumova/velo = \dot{λ}
  variables.nodal_stimulation_current = 40.                           # not used
  variables.vm_value_stimulated = 40.                                 # to which value of Vm the stimulated node should be set (option "valueForStimulatedPoint" of FastMonodomainSolver)
  
elif "Aliev_Panfilov_Razumova_Titin" in variables.cellml_file:   # this is (4, "Titin") in OpenCMISS
  # parameters: I_stim, fiber stretch λ, fiber contraction velocity \dot{λ}
  variables.muscle1_mappings = {
    ("parameter", 0):           ("constant", "Aliev_Panfilov/I_HH"),  # parameter 0 is constant 0 = I_stim
    ("parameter", 1):           ("constant", "Razumova/l_hs"),        # parameter 1 is constant 11 = fiber stretch λ
    ("parameter", 2):           ("constant", "Razumova/rel_velo"),    # parameter 2 is constant 12 = fiber contraction velocity \dot{λ}
    ("connectorSlot", 0): ("state", "Aliev_Panfilov/V_m"),      # expose state 0 = Vm to the operator splitting
    ("connectorSlot", 1): ("algebraic", "Razumova/ActiveStress"),   # expose algebraic 4 = γ to the operator splitting
    ("connectorSlot", 2): ("algebraic", "Razumova/Activation"),     # expose algebraic 5 = α to the operator splitting
  }
  variables.parameters_initial_values = [0, 1, 0]                     # Aliev_Panfilov/I_HH = I_stim, Razumova/l_hs = λ, Razumova/rel_velo = \dot{λ}
  variables.nodal_stimulation_current = 40.                           # not used
  variables.vm_value_stimulated = 40.                                 # to which value of Vm the stimulated node should be set (option "valueForStimulatedPoint" of FastMonodomainSolver)
  
elif "hodgkin_huxley-razumova" in variables.cellml_file:   # this is (4, "Titin") in OpenCMISS
  # parameters: I_stim, fiber stretch λ, fiber contraction velocity \dot{λ}
  variables.muscle1_mappings = {
    ("parameter", 0):           "membrane/i_Stim",          # parameter 0 is I_stim
    ("parameter", 1):           "Razumova/l_hs",            # parameter 1 is fiber stretch λ
    ("connectorSlot", "m1vm"):  "membrane/V",               # expose Vm to the operator splitting
    ("connectorSlot", "m1gout"):"Razumova/activestress",
    ("connectorSlot", "m1alp"): "Razumova/activation",      # expose activation .
    ("connectorSlot", "m1lda"): "Razumova/l_hs",            # fiber stretch λ
  }
  variables.parameters_initial_values = [0, 1]
  variables.nodal_stimulation_current = 40.                           # not used
  variables.vm_value_stimulated = 20.                                 # to which value of Vm the stimulated node should be set (option "valueForStimulatedPoint" of FastMonodomainSolver)

else:
  print("\033[0;31mCellML file {} has no mappings implemented in helper.py\033[0m".format(variables.cellml_file))
  quit()





# load MU distribution and firing times
variables.fiber_distribution = np.genfromtxt(variables.fiber_distribution_file, delimiter=" ", dtype=int)
variables.firing_times = np.genfromtxt(variables.firing_times_file)

# # MU indices in file are 1-based.
# needed_motorneurons = np.max(variables.fiber_distribution[:variables.n_fibers_total])
# # we cannot reindex the MUs so unused MUs are not computed because things like cortical input depend on the MU numbering
# if variables.n_motoneurons < needed_motorneurons:
#     raise Exception(f"Not enough motoneurons: {variables.n_fibers_total} fibers with MUs in [0, {needed_motorneurons}], but n_motoneurons={variables.n_motoneurons}")


# ---------------------
# callback functions
def get_motor_unit_no(fiber_no):
  return int(variables.fiber_distribution[fiber_no % len(variables.fiber_distribution)]-1)

def get_diffusion_prefactor(fiber_no, mu_no):
  return variables.get_conductivity(fiber_no, mu_no) / (variables.get_am(fiber_no, mu_no) * variables.get_cm(fiber_no, mu_no))



