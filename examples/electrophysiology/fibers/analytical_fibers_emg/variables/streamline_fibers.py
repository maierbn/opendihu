
# scenario name for log file
scenario_name = "streamline_fibers"
import sys
import argparse
parser = argparse.ArgumentParser(description='analytical_fibers_emg')

# parameters that will be set by command line arguments
# ------------------------------------------------------
belly_width = 1.0       # this can be set by the command line argument --belly_width, it is available in settings_analytical_fibers_emg_custom_geometry.py as variables.belly_widt
# i.e., run:
# ./analytical_fibers_emg ../settings_analytical_fibers_emg_custom_geometry.py streamline_fibers.py --belly_width=0.1

parser.add_argument('--belly_width', help='Parameter of the geometry', type=float, default=belly_width)

# parse the command line arguments
args,_ = parser.parse_known_args(args=sys.argv[:-2])
locals().update(vars(args))

print("belly_width: {}".format(belly_width))

# timing parameters
# -----------------
end_time = 100.0                    # [ms] end time of the simulation
stimulation_frequency = 100*1e-3    # [ms^-1] sampling frequency of stimuli in firing_times_file, in stimulations per ms, number before 1e-3 factor is in Hertz.
dt_3D = 1                           # [ms] time step width of coupling, when 3D should be performed, also sampling time of monopolar EMG
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

# belly shape function of the muscle, z in [0,1]
def belly(z):
  return 0.5 + (np.sin(z * 2*np.pi - 0.5*np.pi)+1)/2 * belly_width

def map_to_circle(i, j, n_grid_points_x, n_grid_points_y, radius_factor):
  """ 
  map a unit square to a unit circle
  :param i: x coordinate in the unit square, integer, 0 <= i <= n_grid_points_x
  :param j: y coordinate in the unit square, integer, 0 <= j <= n_grid_points_y
  :param n_grid_points_x: number of points in x direction of the unit square
  :param n_grid_points_y: number of points in y direction of the unit square
  :return: (x,y) new point coordinates
  """  
  x = i / (n_grid_points_x-1)
  y = j / (n_grid_points_y-1)
  phi = j / n_grid_points_y * 2.*np.pi
    
  # get segment
  if (n_grid_points_y%2 == 1 and j != int(n_grid_points_y/2.)) or n_grid_points_y%2 == 0:
    if j < n_grid_points_y/2.:   # bottom or side
      if i >= j and i <= n_grid_points_x-1-j:   # bottom
    
        # get layer
        a = ((n_grid_points_y-1)/2.-j)*2.0/(n_grid_points_y-1)    # fraction (in [0,1]) of distance between center of square and line of points
        alpha = float(i - j) / (n_grid_points_x-1-2*j)
        phi = -np.pi/4. + np.pi/2.*alpha
        
        x = np.sin(phi)*a
        y = (-1./np.sqrt(2.) + (-np.cos(phi) + 1./np.sqrt(2.))*a)*a
    
    else:   # top or side
      if i >= n_grid_points_x-1-j and i <= j:   # top
    
        # get layer
        a = (j-(n_grid_points_y-1)/2.)*2.0/(n_grid_points_y-1)    # distance between center of square and line of points
        alpha = float(i - (n_grid_points_y-1-j)) / (2*j-n_grid_points_y+1)
        phi = -np.pi/4 + np.pi/2.*alpha
        
        x = np.sin(phi)*a
        y = (1./np.sqrt(2.) + (np.cos(phi) - 1./np.sqrt(2.))*a)*a
    
  if (n_grid_points_x%2 == 1 and i != int(n_grid_points_x/2.)) or n_grid_points_y%2 == 0:
    if i < n_grid_points_x/2.:   # left
      if j >= i and j <= n_grid_points_y-1-i:   # left
    
        # get layer
        a = ((n_grid_points_x-1)/2.-i)*2.0/(n_grid_points_x-1)    # distance between center of square and line of points
        alpha = float(j - i) / (n_grid_points_y-1-2*i)
        phi = -np.pi/4 + np.pi/2.*alpha
        
        y = np.sin(phi)*a
        x = (-1./np.sqrt(2.) + (-np.cos(phi) + 1./np.sqrt(2.))*a)*a
    
    else:   # right
      if j >= n_grid_points_y-i and j <= i:   # right
    
        # get layer
        a = (i-(n_grid_points_x-1)/2.)*2.0/(n_grid_points_x-1)    # distance between center of square and line of points
        alpha = float(j - (n_grid_points_y-1-i)) / (2*i-n_grid_points_x+1)
        phi = -np.pi/4 + np.pi/2.*alpha
        
        y = np.sin(phi)*a
        x = (1./np.sqrt(2.) + (np.cos(phi) - 1./np.sqrt(2.))*a)*a
     
  # rotate by pi*3/4 
  phi = np.arctan2(y, x) + np.pi/4 + np.pi/2
  r = np.sqrt(x*x + y*y) * radius_factor
  x = np.cos(phi)*r
  y = np.sin(phi)*r
  
  if n_grid_points_x%2 == 1 and i == int(n_grid_points_x/2) and j == int(n_grid_points_y/2):   # center point
    x = 0.
    y = 0.
    
  return (x,y)
  
# loop over points of a structured grid, x,y,z coordinates are i,j,k with size n_nodes_x * n_nodes_y * n_nodes_z
for k in range(n_nodes_z):
  for j in range(n_nodes_y):
    for i in range(n_nodes_x):
      fiber_no = j*n_nodes_x + i
          
      # get belly shape of the muscle
      radius = belly(k/n_nodes_z)
      
      # determine new position
      (x,y) = map_to_circle(i,j,n_nodes_x,n_nodes_y,radius)
      
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
# number of fibers grid
n_custom_fibers_x = 6
n_custom_fibers_y = 6
distance_x = 0.5      # distance between fibers
distance_y = 0.5

n_nodes_per_mm = 1
n_nodes_per_fiber = 100   # for a 10cm muscle 100 nodes
n_custom_fibers = 0     # counter for actually generated fibers

# loop over fibers in a 2D grid
for j in range(n_custom_fibers_y):
  for i in range(n_custom_fibers_x):
    
    x_center = i*distance_x - (n_custom_fibers_y-1)/2*distance_x
    y_center = j*distance_y - (n_custom_fibers_x-1)/2*distance_y
    
    # loop over node positions of the current fiber
    node_positions = []
    for k in range(0,n_nodes_per_fiber):
      
      # get belly shape of the muscle
      radius = belly(k/n_nodes_per_fiber)
      radius *= 0.8                     # no fibers at the outer 20% of the muscle
      
      # determine position
      (x,y) = map_to_circle(i, j, n_custom_fibers_x, n_custom_fibers_y, radius)
      z = k/n_nodes_per_fiber * 10      # 10 cm in z direction
      
      # add position to list of node positions for current fiber
      node_positions.append([x,y,z])
    
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
enable_matrix_output = True

firing_times_file = "../../../input/MU_firing_times_often.txt"                # this file specifies when the motor units fire
fiber_distribution_file = "../../../input/MU_fibre_distribution_10MUs.txt"    # this file specifies which fibers are in which motor units
