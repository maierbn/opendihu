# parallel fiber estimation, Laplace 3D
#

import numpy as np
import sys, os
import argparse
import distutils.util
sys.path.append(os.path.dirname(__file__))
import compute_sizes

# parse arguments
rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

# define command line arguments
mbool = lambda x:bool(distutils.util.strtobool(x))   # function to parse bool arguments
parser = argparse.ArgumentParser(description='Program that generates 3D meshes and 1D fiber meshes out of the surface of a muscle. The program implements the parallel algorithm involving partitioning, harmonic maps and streamline tracing.')
parser.add_argument('--input_filename_or_splines_or_stl',             default="splines",  help='The input filename of the surface triangulation. "splines" and "stl" will be replaced by existing files in input directory')
parser.add_argument('--output_filename',                              default="",         help='Filename of the output bin file, "0x0" gets replaced by the actual number of fibers in the file.')
parser.add_argument('--bottom_z_clip',                    type=float, default=72.0,       help='bottom z value of the muscle volume where potential flow will be computed')
parser.add_argument('--top_z_clip',                       type=float, default=220.0,      help='top z value of the muscle volume')
parser.add_argument('--element_size',                     type=float, default=0.1,        help='Size of one element of a fiber, i.e. the distance between the points on the fiber in the result')
parser.add_argument('--refinement_factor',                type=int,   default=2,          help='Factor, used in x,y,z direction by which the mesh should be refined prior to solving the laplace problem and tracing the streamlines. The number of elements is increased by 8^refinement_factor')
parser.add_argument('--improve_mesh',                     type=mbool, default=True,       help='If Laplacian smoothing and fixing of invalid quadrilaterals should be enabled. This increases the runtime a bit.')
parser.add_argument('--use_gradient_field',               type=mbool, default=False,      help='If the gradient field should be computed explicitely and used for tracing the fibers. If false, only solution values are computed and the gradient directions are computed from the gradient values of the ansatz functions, which is a bad idea for linear ansatz functions.')
parser.add_argument('--use_neumann_bc',                   type=mbool, default=False,      help='True = use Neumann BC, False = use Dirichlet BC for the potential flow Laplace problem')
parser.add_argument('--n_elements_z_per_subdomain', '-z', type=int,   default=50,         help='Number of elements in z direction per subdomain.')
parser.add_argument('--n_elements_x_per_subdomain', '-x', type=int,   default=4,          help='Number of elements in x direction per subdomain.')
parser.add_argument('--n_fine_grid_fibers', '-m',         type=int,   default=0,          help='Number of fine grid fibers to interpolate between the key fibers, parameter is called m.')
parser.add_argument('--max_level', '-l',                  type=int,   default=4,          help='Maximum recursion level l, the required number of processes is 8^l. The maximum recursion level is also reached when there are not enough processes left, so a too high value for l is no problem.')
parser.add_argument('--ghost_layer_width',                type=int,   default=-1,         help='Thickness of the ghost layer in number of elements. The value -1 means equal to 4*refinement_factor.')
parser.add_argument('--max_area_factor',                  type=float, default=100.0,      help='factor only for triangulation_type 1, approximately the minimum number of triangles that will be created because of a maximum triangle area constraint.')
parser.add_argument('--program_name',                                 default="generate", help='If the program is run with the quadratic ansatz functions, this information is only for the scenario name.')

# parse command line arguments and assign values to variables module
args = vars(parser.parse_args(args=sys.argv[:-2]))

# set input_mesh_name
use_splines = False
if args["input_filename_or_splines_or_stl"] == "splines":
  input_mesh_filename = "../../../electrophysiology/input/biceps.surface.pickle"
  use_splines = True
elif args["input_filename_or_splines_or_stl"] == "stl":
  input_mesh_filename = "../../../electrophysiology/input/biceps_splines.stl"
else:
  input_mesh_filename = args["input_filename_or_splines_or_stl"]

bottom_z_clip = args["bottom_z_clip"]
top_z_clip = args["top_z_clip"]
element_size = args["element_size"]
refinement = args["refinement_factor"]
improve_mesh = args["improve_mesh"]
use_gradient_field = args["use_gradient_field"]
use_neumann_bc = args["use_neumann_bc"]
use_quadratic = "quadratic" in args["program_name"]
ghost_layer_width = args["ghost_layer_width"]
if ghost_layer_width == -1:
  ghost_layer_width = 2*refinement
n_elements_z_per_subdomain = args["n_elements_z_per_subdomain"]
n_elements_x_per_subdomain = args["n_elements_x_per_subdomain"]
n_fine_grid_fibers = args["n_fine_grid_fibers"]
max_area_factor = args["max_area_factor"]
max_level = args["max_level"]
new_max_level = min(max_level,(int)(np.log(n_ranks) / np.log(8)))
if new_max_level != max_level:
  if rank_no == 0:
    print("Adjusting max_level (l) from {} to {} because there are only {} processes and 8^{} = {}, 8^{} = {}".
      format(max_level, new_max_level, n_ranks, max_level, 8**max_level, new_max_level, 8**new_max_level))
  max_level = new_max_level

scenario_name = "l{}_m{}_n{}_{}{}{}{}{}_{}".format(max_level, n_fine_grid_fibers, n_elements_x_per_subdomain, 
  "q" if use_quadratic else "l", 
  "N" if use_neumann_bc else "D", 
  refinement, 
  "g" if use_gradient_field else "s", 
  "splines" if use_splines else "stl",
  "i" if improve_mesh else "n")

if args["output_filename"] == "":
  output_filename = "0x0fibers_{}.bin".format(scenario_name)
else:
  output_filename = args["output_filename"]

n_nodes_per_fiber = (top_z_clip - bottom_z_clip) / element_size
n_nodes_per_fiber = 2 * (n_nodes_per_fiber // 2) + 1   # make number odd

if rank_no == 0:
  print("Parameters:")
  print("  splines_or_stl:      {}".format(args["input_filename_or_splines_or_stl"]))
  print("  input_mesh_filename: {}".format(input_mesh_filename))
  print("  output_filename:     {}".format(output_filename))
  print("  bottom_z_clip:       {}".format(bottom_z_clip))
  print("  top_z_clip:          {}".format(top_z_clip))
  print("  element_size:        {}".format(element_size))
  print("  n_nodes_per_fiber:   {}".format(n_nodes_per_fiber))
  print("  refinement:          {}".format(refinement))
  print("  improve_mesh:        {}".format(improve_mesh))
  print("  use_gradient_field:  {}".format(use_gradient_field))
  print("  use_neumann_bc:      {}".format(use_neumann_bc))
  print("  use_quadratic:       {} (program name: \"{}\")".format(use_quadratic,args["program_name"]))
  print("  ghost_layer_width:   {}".format(ghost_layer_width))
  print("  scenario_name:       {}".format(scenario_name))
  print("  n_elements_x_per_subdomain:  {}".format(n_elements_x_per_subdomain))
  print("  n_elements_z_per_subdomain:  {}".format(n_elements_z_per_subdomain))
  print("  n_fine_grid_fibers:  m = {}".format(n_fine_grid_fibers))
  print("  max_level:           l = {}".format(max_level))
  compute_sizes.output_size(n_fine_grid_fibers, max_level, n_elements_x_per_subdomain)

config = {
  "scenarioName": scenario_name,
  "solverStructureDiagramFile": None,
  "mappingsBetweenMeshesLogFile": None,
  "logFormat": "csv",
  "solverStructureDiagramFile": "solver_structure.txt",
  "mappingsBetweenMeshesLogFile": "",
  "Solvers": {
    "linearSolver": {
      "relativeTolerance":  1e-4,
      "absoluteTolerance":  1e-10,         # 1e-10 absolute tolerance of the residual    
      "maxIterations":      5e3,
      "solverType":         "gmres",
      "preconditionerType": "sor",
      "dumpFormat":         None,
      "dumpFilename":       None,
    }
  },
  "ParallelFiberEstimation" : {
    "inputMeshFilename":          input_mesh_filename,   # this is the input filename
    "resultFilename":             output_filename,       # this is the output filename, the numbers <a>x<b> are adjusted automatically
    
    "bottomZClip":                bottom_z_clip,         # 82 (72), bottom z value of the muscle volume to simulate the potential flow in
    "topZClip":                   top_z_clip,            # 250 (220), top z value of the muscle volume
    "finalBottomZClip":           bottom_z_clip,         # 82 (72), bottom z value of the final fibers, fibers will be cropped and resampled to nNodesPerFiber between finalBottomZClip and finalTopZClip
    "finalTopZClip":              top_z_clip,            # 250 (220), top z value of the final fibers, fibers will be cropped and resampled to nNodesPerFiber between finalBottomZClip and finalTopZClip
    "useNeumannBoundaryConditions": use_neumann_bc,             # which type of boundary conditions at top and bottom should be used, Neumann or Dirichlet type  
    "nElementsXPerSubdomain":     n_elements_x_per_subdomain,   # 4 number of elements in x and y-direction per subdomain
    "nElementsZPerSubdomain":     n_elements_z_per_subdomain,   # 50 number of elements in z-direction per subdomain
    "nFineGridFibers":            n_fine_grid_fibers,           # number of additional fine fibers that are interpolated between the main "key" fibers, the key fibers are traced
    "useGradientField":           use_gradient_field,           # set to False
    "maxLevel":                   max_level,                    # maximum level (0=1 process, 1=8 processes, 2=64 processes)
    "lineStepWidth":              0.01,                  # line width for tracing of fibers
    "nNodesPerFiber": n_nodes_per_fiber,                 # number of nodes in each final fiber
    "maxIterations":              1e5,                   # maximum number of iterations in the solver
    "maxAreaFactor":              max_area_factor,       # factor only for triangulation_type 1, approximately the minimum number of triangles that will be created because of a maximum triangle area constraint
    
    "improveMesh":                improve_mesh,          # smooth the 2D meshes, required for bigger meshes or larger amount of ranks
    "refinementFactors": [refinement,refinement,refinement],         # [2,2,2] factors in x,y,z direction by which the mesh should be refined prior to solving the laplace problem and tracing the streamlines
    "laplacianSmoothingNIterations": 10,                 # number of Laplacian smoothing iterations on the final fibers grid
    "ghostLayerWidth":            ghost_layer_width,     # width of the ghost layer of elements that is communicated between the subdomains, such that boundary streamlines do not leave the subdomains during tracing
    
    "FiniteElementMethod": {
      "meshName":   "potentialFlow",
      "solverName": "linearSolver",
      "dirichletBoundaryConditions": {},    # the boundary conditions will be set automatically by the program
      "dirichletOutputFilename":  None,     # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
      "neumannBoundaryConditions":   [],
      "prefactor": 1.0,
      "maxIterations": 1e5,
      "inputMeshIsGlobal": True,
      "slotName": "",
    },
    "OutputWriter": [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/biceps", "binary": True, "fixedFormat": False, "combineFiles": False, "fileNumbering": "incremental"},
    ]
  }
}
