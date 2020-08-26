# Linear elasticity with active stress and artifical electrophysiology on fibers.
# Use ParaView for visualization. Display fibers as "Point Gaussian".
# Add a Glyph filter with arrows for field "-rhsNeumannBC", this is the applied traction.
#
# This settings file works for linear or quadratic ansatz functions, 
# if the computation of number of nodes, mx,my,mz (lines ~30) are changed accordingly.

import numpy as np
import sys, os

# number of 3D elements in x, y and z direction, z is the direction of contraction
nx = 3
ny = 3
nz = 5      # direction of muscle
physical_extent = [3, 3, 10]    # extent of the 3D mesh [x,y,z]

# timestep width
dt = 0.1

# number of fibers in x and y direction
n_fibers_x = 5
n_fibers_y = 5
n_nodes_per_fiber = 20        # number of nodes per fiber


# number of nodes of the 3D element
# for quadratic ansatz functions
#mx = 2*nx+1
#my = 2*ny+1
#mz = 2*nz+1

# for linear ansatz functions
mx = nx+1
my = ny+1
mz = nz+1

# get the number of MPI ranks and the own rank no. These two arguments are passed to this script by the c++ program as the last two arguments.
rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

def set_artifical_activation_values(n_dofs_global, n_nodes_global_per_coordinate_direction, time_step_no, current_time, values, global_natural_dofs, fiber_index):
  # n_dofs_global:       (int) global number of dofs in the mesh where to set the values
  # n_nodes_global_per_coordinate_direction (list of ints)   [mx, my, mz] number of global nodes in each coordinate direction. 
  #                       For composite meshes, the values are only for the first submesh, for other meshes sum(...) equals n_dofs_global
  # time_step_no:        (int)   current time step number
  # current_time:        (float) the current simulation time
  # values:              (list of floats) all current local values of the field variable, if there are multiple components, they are stored in struct-of-array memory layout 
  #                       i.e. [point0_component0, point0_component1, ... pointN_component0, point1_component0, point1_component1, ...]
  # global_natural_dofs  (list of ints) for every local dof no. the dof no. in global natural ordering
  # additional_argument: The value of the option "additionalArgument", can be any Python object.

  # loop over local dofs
  for local_dof_no in range(len(values)):
    global_dof_no = global_natural_dofs[local_dof_no]
    
    # node index along fiber, k ∈ [0,N)
    k = global_dof_no
    N = n_nodes_global_per_coordinate_direction[0]    # number of points per fiber
    
    # fiber index (x,y) ∈ [0,n_fibers_x) x [0,n_fibers_y)
    fiber_x = fiber_index % n_fibers_x
    fiber_y = int(fiber_index / n_fibers_x)
    
    # current time, t ∈ [0,T)
    t = current_time
    T = 10.0
    
    # set activation value ∈ [0,1]
    values[local_dof_no] = np.sin((0.25*k/N + 0.1*fiber_x/n_fibers_x + 0.1*fiber_y/n_fibers_y + t/T) * 2*np.pi) ** 2


# ---------------
# define boundary conditions of the elasticity problem (for linear elements)
# fix muscle at the top (z+)
dirichlet_bc = {}

# top plane (z+)
for j in range(0,my):
  for i in range(0,mx):
    k = mz-1
    dirichlet_bc[k*mx*my + j*mx + i] = [None, None, 0.0]      # [x,y,z], fix only dof in z direction

# front edge on top plan (y-)
for i in range(0,mx):
  k = mz-1
  j = 0
  dirichlet_bc[k*mx*my + j*mx + i][1] = 0     # fix dof in y direction

# left edge on top plane (x-)
for j in range(0,my):
  k = mz-1
  i = 0
  dirichlet_bc[k*mx*my + j*mx + i][0] = 0     # fix dof in x direction

# traction boundary conditions, pull downwards (negative z direction)
k = nz-1
neumann_bc = [{"element": k*nx*ny + j*nx + i, "constantVector": [0.0,0.0,-3.0], "face": "2-"} for j in range(ny) for i in range(nx)]

#dirichlet_bc = {}
#neumann_bc = []

# ---------------
# create meshes
meshes = {
  "elasticityMesh": {
    "nElements":          [nx, ny, nz],             # number of elements in the 3 coordinate directions  
    "inputMeshIsGlobal":  True,                     # if the specifiation of the mesh is given in global numbers
    "physicalExtent":     physical_extent,          # the physical dimensions of the mesh, in this case the same as the number of elements
    "physicalOffset":     [0, 0, 0],                # lower left bottom coordinates of the mesh
  }
}

# add fiber meshes

# loop over fibers
for j in range(n_fibers_y):
  for i in range(n_fibers_x):
    fiber_index = j*n_fibers_x + i
  
    # determine node positions of the fiber mesh
    node_positions = []
    start_point = np.array([(i+1) / (n_fibers_x+1) * physical_extent[0], (j+1) / (n_fibers_y+1) * physical_extent[1], 0])
    hz = physical_extent[2] / (n_nodes_per_fiber-1)
    
    for k in range(n_nodes_per_fiber):
      node_positions.append(list(start_point + np.array([0,0,k*hz])))
  
    # add fiber mesh
    name = "MeshFiber_{}".format(fiber_index)
    meshes[name] =  {
      "nElements":          n_nodes_per_fiber-1,      # number of elements of the fiber 
      "inputMeshIsGlobal":  True,                     # if the specifiation of the mesh is given in global numbers
      "nodePositions":      node_positions,           # node positions of this fiber mesh
      "nodeDimension":      3,                        # only needed if the nodePositions list is a consecutive array [x0,y0,z0,x1,y1,z1,...], but here it is not needed as we have the proper nodes [[x0,y0,z0], [x1,y1,z1], [...], ...]
    }
    
# print meshes for debugging
if False:
  import pprint
  pp = pprint.PrettyPrinter(indent=4)
  pp.pprint(meshes)
    
config = {
  "scenarioName":                   "linear_elasticity",             # scenario name for the log file
  "logFormat":                      "csv",            # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     "solver_structure.txt",          # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile":   "mappings_between_meshes.txt",   # log file for mappings between meshes
  "Meshes":                 meshes,
  "MappingsBetweenMeshes":  {"MeshFiber_{}".format(i) : "elasticityMesh" for i in range(n_fibers_x*n_fibers_y)},    # initialize mapping between fiber meshes and 3D mesh
  "Solvers": {
    "linearElasticitySolver": {   # solver for linear elasticity
      "relativeTolerance":  1e-15,                    # relative tolerance for the solver, this relates to the current residual norm divided by the norm of the rhs
      "absoluteTolerance":  1e-10,                    # absolute tolerance of the residual        
      "solverType":         "gmres",                  # which PETSc linear system solver to use
      "preconditionerType": "none",                   # which PETSc preconditioner to use
      "maxIterations":      1e4,                      # maximum number of iterations of the linear solver
      "dumpFormat":         "matlab",                 # "default", "ascii", "matlab", output format type for system matrix and right hand side
      "dumpFilename":       "",                       # filename for output of system matrix and right hand side after every solve, set to "" to disable
    }, 
  },
  "Coupling": {
    "timeStepWidth":          dt,                     # time step width of the coupling, i.e. when the diffusion and reaction terms will be evaluated
    "logTimeStepWidthAsKey":  "dt_coupling",          # the key under which the timeStepWidth value will be written to the log file
    "durationLogKey":         "duration_total",       # a key to find the total computational duration in the log file
    "timeStepOutputInterval": 10,                     # how often to print information about iteration to the console
    "endTime":                10,                     # end time of the simulation
    
    "connectedSlotsTerm1To2": {0:0},                  # connect output slots of the two coupled terms: artifical_activation (1D fibers) -> activation (3D elasticity) (also see the solver_structure.txt file)
    "connectedSlotsTerm2To1": None,                   # no connection from 3D to 1D
 
    # artifical electrophysiology
    "Term1": {  
      "MultipleInstances": {
        "logKey":                     "duration_fibers",        # a key for the log file to find the duration of fibers
        "ranksAllComputedInstances":  list(range(n_ranks)),     # list of all MPI ranks that do this computation
        "nInstances":                 n_fibers_x * n_fibers_y,  # number of instances (fibers)
        "instances": 
        [{
          "ranks": list(range(n_ranks)),              # list of MPI ranks that compute this particular fiber
              
          # dummy solver that sets the activation field variable by a python callback
          "PrescribedValues": {
            "meshName":               "MeshFiber_{}".format(fiber_index),    # the mesh to use for the field variables
            "numberTimeSteps":        1,              # number of timesteps to call the callback functions subsequently, this is usually 1 for prescribed values, because it is enough to set this term only once per time step
            "timeStepOutputInterval": 20,             # if the time step should be written to console, a value > 10 produces no output
            "slotNames":              [],             # names of the connector slots, if the global option "connectedSlots" is used
            
            # a list of field variables that will get values assigned in every timestep, by the provided callback function
            "fieldVariables1": [
              {"name": "artifical_activation", "callback": set_artifical_activation_values},
            ],
            "fieldVariables2":     [],
            "additionalArgument":  fiber_index,       # a custom argument to the fieldVariables callback functions, this will be passed on as the last argument. Here it is the fiber index
          }   
        } for fiber_index in range(n_fibers_x * n_fibers_y)],
    
        # output of elasticity values is disabled, because they are proportional to the active stress which is output by the linear elasticity solver
        "OutputWriter" : [
          {"format": "Paraview", "outputInterval": 1, "filename": "out/fibers", "binary": True, "fixedFormat": True, "onlyNodalValues":True, "combineFiles": True, "fileNumbering": "incremental"}
        ]
      }
    },
    
    # linear elasticity
    "Term2": {      
    
      # timestepping solver that solves quasi-static elasticity formulation
      "QuasiStaticLinearElasticitySolver": {
        "fiberDirection":           [0, 0, 1],      # direction for anisotropy of elasticity formulation
        "PotentialFlow":            None,           # if fiberDirection is not set to a constant direction, a potential flow simulation can be used where the fiber direction is set to the streamlines of the flow through the volume. In this case, set "PotentialFlow" to the settings for the FEM for the potential flow.
        "maximumActiveStress":      1.0,            # scaling factor to the active stress, σ_active = activation * anisotropyTensor * maximumActiveStress
        "strainScalingCurveWidth":  1.0,            # parameter for strain-stress curve of active stress, has no effect, because strain-stress curve is commented out in the code
        "scalingFactor":            1.0,            # scaling factor for displacements, to overrate them, if != 0 it is only for visualization purposes and not physical
        "inputMeshIsGlobal":        True,           # if boundary conditions are specified in global numbering
        "slotNames":                [],             # names of the connector slots, if the global option "connectedSlots" is used
          
        # Anisotropy for active stress.
        # The tensor is given in a local basis where the fiber direction is (1,0,0), one list item = same tensor for all elements, multiple list items = a different tensor for each element.
        # The tensor has to be symmetric.
        "anisotropyTensor": [              
          1, 0, 0,
          0, 0, 0,
          0, 0, 0
        ],
        
        # linear elasticity finite element method
        "FiniteElementMethod" : {   
          "meshName":             "elasticityMesh",       # reference to the mesh, defined under "Meshes"
          "solverName":           "linearElasticitySolver",   # reference to the linear system solver, defined under "Solvers"
          "prefactor":            1.0,                    # prefactor c of the equation c*Δu = f
          "inputMeshIsGlobal":    True,                   # if boundary conditions are specified in global numbering
          "dirichletBoundaryConditions":  dirichlet_bc,   # dirichlet boundary conditions
          "dirichletOutputFilename":      None,           # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
      
          "neumannBoundaryConditions":    neumann_bc,     # neumann boundary conditions
          "divideNeumannBoundaryConditionValuesByTotalArea": False,  # if the neumann boundary condition vectors should be divided by the total surface area where surface loads are applied, this allows to specify the total force that acts on the surface. If set to False (default), the given traction is a per-surface quantity.
          "slotName":             "",                     # name of the connector slot of the solution field variable, needed if the global option "connectedSlots" is used
                
          # material parameters
          "bulkModulus":        150,                      # bulk modulus, K, material parameter for compressibility
          "shearModulus":       20,                       # shear modulus, μ
        },
        
        "OutputWriter": [
          {"format": "Paraview", "outputInterval": 1, "filename": "out/elasticity", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
          #{"format": "PythonFile", "filename": "out/elasticity", "outputInterval": 1, "binary":False, "onlyNodalValues":True, "fileNumbering": "incremental"},
        ]
      }
    }
  }
}

# print config for debugging
if False:
  import pprint
  pp = pprint.PrettyPrinter(indent=2)
  pp.pprint(config["Coupling"])
