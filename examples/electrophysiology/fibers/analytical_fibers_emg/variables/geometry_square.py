
# scenario name for log file
scenario_name = "custom_geometry"

# timing parameters
# -----------------
end_time = 100.0                    # [ms] end time of the simulation
stimulation_frequency = 100*1e-3    # [ms^-1] sampling frequency of stimuli in firing_times_file, in stimulations per ms, number before 1e-3 factor is in Hertz.
dt_3D = 1e-1                           # [ms] time step width of coupling, when 3D should be performed, also sampling time of monopolar EMG
output_timestep = 1                 # [ms] timestep for output surface EMG, 0.5
output_timestep_fibers = 1          # [ms] timestep for fiber output, 0.5
output_timestep_big = 1             # [ms] timestep for output files of 3D intramuscular EMG data

# custom 3D mesh
# ----------------------
import numpy as np
custom_meshes = {}
node_positions = []

# number of nodes of the 3D mesh
n_nodes_x = 6
n_nodes_y = 6
n_nodes_z = 100

# belly shape function of the muscle
def belly(r,z):
  return r*(0.5 + (np.sin(z * 2*np.pi - 0.5*np.pi)+1)/2)

# loop over points of a structured grid, x,y,z coordinates are i,j,k with size n_nodes_x * n_nodes_y * n_nodes_z
for k in range(n_nodes_z):
  for j in range(n_nodes_y):
    for i in range(n_nodes_x):
      fiber_no = j*n_nodes_x + i
        
      x = i-(n_nodes_x-1)/2
      y = j-(n_nodes_y-1)/2
      
      # get polar coordinates of current point
      phi = np.arctan2(y,x)
      r = np.linalg.norm(np.array((x,y)))
          
      # get belly shape of the muscle
      r_new = belly(r,k/n_nodes_z)
      
      # determine new position
      x = r_new*np.cos(phi)
      y = r_new*np.sin(phi)
      z = k/n_nodes_z * 10      # 10 cm in z direction
        
      node_positions.append([x,y,z])
    
print("n nodePositions: {}".format(len(node_positions)))
print("nElements: {}".format((n_nodes_x-1) * (n_nodes_y-1) * (n_nodes_z-1)))
    
custom_meshes["3Dmesh"] = {
  "nElements": [(n_nodes_x-1), (n_nodes_y-1), (n_nodes_z-1)],
  "nodePositions": node_positions,
  "inputMeshIsGlobal": True,
  "setHermiteDerivatives": False,
  "logKey": "3Dmesh",
  "nRanks": [1,1,1],
}

# custom 1D fiber meshes
# ----------------------
# number of fiber grid, only those fibers that are actually inside the muscle are generated
n_custom_fibers_x = 8
n_custom_fibers_y = 8
n_nodes_per_fiber = 1000
alpha = 0.4*np.pi       # angles to skew the fibers
beta = 0.2*np.pi

n_custom_fibers = 0     # counter for actually generated fibers

# loop over fibers in a 2D grid
for j in range(n_custom_fibers_y):
  for i in range(n_custom_fibers_x):
    
    x_center = i-(n_custom_fibers_y-1)/2
    y_center = j-(n_custom_fibers_x-1)/2
    
    # loop over node positions of the current fiber
    node_positions = []
    for k in range(-n_nodes_per_fiber//2,n_nodes_per_fiber//2):
      
      # determine position
      x = x_center + k/n_nodes_per_fiber * np.tan(alpha)
      y = y_center + k/n_nodes_per_fiber * np.tan(beta)
      z = (k/n_nodes_per_fiber + 0.5) * 10      # 10 cm in z direction
      
      # determine if fiber point is inside muscle volume
      r_base = np.linalg.norm(np.array((i-(n_nodes_x-1)/2,j-(n_nodes_y-1)/2)))
      r = np.linalg.norm(np.array((x,y)))
      
      # if point is outside muscle, by using the belly heuristic
      if r > belly(r_base,k/n_nodes_per_fiber + 0.5):        
        if len(node_positions) == 0:            # if fiber inside the muscle has not yet started, continue, waiting for the first point inside the muscle
          continue
        else:               # if there were already points inside the muscle, finish the current fiber
          break
      
      # add position to list of node positions for current fiber
      node_positions.append([x,y,z])
    
    # if there were node positions inside the muscle, add fiber to dict of fibers
    if len(node_positions) == 0:
      print("Fiber ({},{}) is completely outside the 3D mesh!".format(i,j))
    else:          
      custom_meshes["MeshFiber_{}".format(n_custom_fibers)] = {
        "nElements": len(node_positions)-1,
        "nodePositions": node_positions,
        "inputMeshIsGlobal": True,
        "setHermiteDerivatives": False,
        "nRanks": [1],
      }
      n_custom_fibers += 1
      
# other options
paraview_output = True
adios_output = False
exfile_output = False
python_output = False
enable_surface_emg = True
disable_firing_output = False

#fiber_file = "../../../input/left_biceps_brachii_13x13fibers.bin"
fiber_file = "../../../input/left_biceps_brachii_7x7fibers.bin"
firing_times_file = "../../../input/MU_firing_times_real.txt"
fiber_distribution_file = "../../../input/MU_fibre_distribution_10MUs.txt"
