# Hyperelastic custom material, static problem
import numpy as np
import sys, os

scenario_name = "static"

# number of elements
nx = 2    # 2
ny = 2    # 2
nz = 5    # 5

# size of box
physical_extent = [2, 2, 5]   # [cm]

# number of nodes (for quadratic elements)
mx = 2*nx + 1
my = 2*ny + 1
mz = 2*nz + 1

# Boundary conditions
# Dirichlet boundary conditions
dirichlet_bc = {}

xpos = 0.0
ypos = 0.0
zpos = 0.0

# parse own MPI rank no and total number of processes
own_rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

# bottom plane
for j in range(0,my):
  for i in range(0,mx):
    dirichlet_bc[j*(mx) + i] = [None,None,zpos]

if False:
  # vertical plane at x=0
  for k in range(0,mz):
    for j in range(0,my):
      dirichlet_bc[k*mx*my + j*mx + 0] = [xpos,None,None]

  # vertical plane at y=0
  for k in range(0,mz):
    for i in range(0,mx):
      dirichlet_bc[k*mx*my + i] = [None,ypos,None]

  # vertical edge at x=0, y=0
  for k in range(0,mz):
    dirichlet_bc[k*mx*my] = [xpos,ypos,None]

if True:
  # edge z=0, y=0
  for i in range(0,mx):
    dirichlet_bc[i] = [None,ypos,zpos]

if True:
  # edge, z=0, x=0
  for j in range(0,my):
    dirichlet_bc[j*(mx)] = [xpos,None,zpos]

# bottom corner
dirichlet_bc[0] = [xpos,ypos,zpos]
dirichlet_bc[1] = [None,ypos,zpos]

# Neumann boundary conditions (traction)
# constant forces at right plane (2+ = positive z direction)
# [N]
neumann_bc = [{"element": (nz-1)*nx*ny + j*nx + i, "constantVector": [0,0,10], "face": "2+"} for j in range(ny) for i in range(nx)]

#dirichlet_bc = {}
#neumann_bc = []

# material parameters
c1 = 5        # [N/cm^2] first Mooney-Rivlin parameter
c2 = 7.5      # [N/cm^2] second Mooney-Rivlin parameter
c =  10       # [-]      bulk modulus (high = less compressible)
material_parameters = [c1, c2, c]

# load
#constant_body_force = (0,0,-9.81e-4)   # [cm/ms^2], gravity constant for the body force in negative z direction
constant_body_force = (0,0,0)
    
config = {
  "scenarioName": scenario_name,            # scenario name to find the run in the log file
  "logFormat":    "csv",                    # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     None,   # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile":   None,   # log file for mappings between meshes
  
  "Meshes": {
    "3Dbox_quadratic": {
      "inputMeshIsGlobal":          True,                       # boundary conditions are specified in global numberings, whereas the mesh is given in local numberings
      "nElements":                  [nx, ny, nz],               # number of quadratic elements in x, y and z direction
      "physicalExtent":             physical_extent,            # physical size of the box
      "physicalOffset":             [0, 0, 0],                  # offset/translation where the whole mesh begins
    }
  },
  
  "Solvers": {
    "mechanicsSolver": {                  
      "relativeTolerance":  1e-5,           # 1e-10 relative tolerance of the linear solver
      "absoluteTolerance":  1e-10,          # 1e-10 absolute tolerance of the residual of the linear solver
      "solverType":         "lu",           # type of the linear solver: cg groppcg pipecg pipecgrr cgne nash stcg gltr richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs fbcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly qcg bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres tsirm cgls
      "preconditionerType": "none",         # type of the preconditioner
      "maxIterations":       1e4,           # maximum number of iterations in the linear solver
      "snesMaxFunctionEvaluations": 1e8,    # maximum number of function iterations
      "snesMaxIterations":   10,            # maximum number of iterations in the nonlinear solver
      "snesRelativeTolerance": 1e-5,        # relative tolerance of the nonlinear solver
      "snesAbsoluteTolerance": 1e-5,        # absolute tolerance of the nonlinear solver
      "snesLineSearchType": "l2",           # type of linesearch, possible values: "bt" "nleqerr" "basic" "l2" "cp" "ncglinear"
      "snesRebuildJacobianFrequency": 2,    # how often the jacobian should be recomputed, -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time the Jacobian is built etc. -2 means rebuild at next chance but then never again 
      "dumpFilename":        "",            # filename to dump system matrix and right hand side after every solve, set to "" to disable
      "dumpFormat":          "matlab",      # format of dump, one of default, ascii, matlab
    }
  },
  
  "HyperelasticitySolver": {
    "durationLogKey":             "nonlinear",               # key to find duration of this solver in the log file
    
    "materialParameters":         material_parameters,       # material parameters of the Mooney-Rivlin material
    "displacementsScalingFactor": 1.0,                       # scaling factor for displacements, only set to sth. other than 1 only to increase visual appearance for very small displacements
    "residualNormLogFilename":    "log_residual_norm.txt",   # log file where residual norm values of the nonlinear solver will be written
    "useAnalyticJacobian":        True,                      # whether to use the analytically computed jacobian matrix in the nonlinear solver (fast)
    "useNumericJacobian":         False,                     # whether to use the numerically computed jacobian matrix in the nonlinear solver (slow), only works with non-nested matrices, if both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian
      
    "dumpDenseMatlabVariables":   False,                     # whether to have extra output of matlab vectors, x,r, jacobian matrix (very slow)
    # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables all all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian
    
    # mesh
    "meshName":                   "3Dbox_quadratic",         # name of the 3D mesh, it is defined under "Meshes" at the beginning of this config
    "fiberMeshNames":             [],                        # fiber directory for anisotropic material, fiber meshes that will be used to determine the fiber direction, for multidomain there are no fibers so this would be empty list
    "fiberDirection":             [0,0,1],                   # fiber directory for anisotropic material, if fiberMeshNames is empty, directly set the constant fiber direction, in element coordinate system
 
    # solver
    "solverName":                 "mechanicsSolver",         # name of the nonlinear solver configuration, it is defined under "Solvers" at the beginning of this config
    #"loadFactors":                [0.25, 0.66, 1.0],        # load factors for every timestep
    "loadFactors":                [],                        # no load factors, solve problem directly
    "loadFactorGiveUpThreshold":   0.25,                     # if the solution diverges, smaller load factors are tried automatically. This option is the threshold when to finally abort the solve instead of trying a smaller load factor.
    "nNonlinearSolveCalls":       1,                         # how often the nonlinear solve should be repeated (should always be one)
    
    # boundary and initial conditions
    "dirichletBoundaryConditions": dirichlet_bc,                        # the initial Dirichlet boundary conditions that define values for displacements u and velocity v
    "neumannBoundaryConditions":   neumann_bc,                          # Neumann boundary conditions that define traction forces on surfaces of elements
    "divideNeumannBoundaryConditionValuesByTotalArea": False,           # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
    "updateDirichletBoundaryConditionsFunction": None,                  # function that updates the dirichlet BCs while the simulation is running
    "updateDirichletBoundaryConditionsFunctionCallInterval": 1,         # every which step the update function should be called, 1 means every time step
    "inputMeshIsGlobal":           True,                                # boundary conditions are given as global numbers
    
    "initialValuesDisplacements":  [[0.0,0.0,0.0] for _ in range(mx*my*mz)],     # the initial values for the displacements, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
    "initialValuesVelocities":     [[0.0,0.0,0.0] for _ in range(mx*my*mz)],     # the initial values for the velocities, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
    "extrapolateInitialGuess":     True,                                # if the initial values for the dynamic nonlinear problem should be computed by extrapolating the previous displacements and velocities
    "constantBodyForce":           constant_body_force,                 # a constant force that acts on the whole body, e.g. for gravity
    
    "dirichletOutputFilename":     "out/"+scenario_name+"/dirichlet_boundary_conditions",          # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    
    # define which file formats should be written
    # 1. main output writer that writes output files using the quadratic elements function space. Writes displacements, velocities and PK2 stresses.
    "OutputWriter" : [
      
      # Paraview files
      {"format": "Paraview", "outputInterval": 1, "filename": "out/"+scenario_name+"/u", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      
      # Python callback function "postprocess"
      #{"format": "PythonCallback", "outputInterval": 1, "callback": postprocess, "onlyNodalValues":True, "filename": ""},
    ],
    # 2. additional output writer that writes also the hydrostatic pressure
    "pressure": {   # output files for pressure function space (linear elements), contains pressure values, as well as displacements and velocities
      "OutputWriter" : [
        {"format": "Paraview", "outputInterval": 1, "filename": "out/"+scenario_name+"/p", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      ]
    },
    # 4. output writer for debugging, outputs files after each load increment, the geometry is not changed but u and v are written
    "LoadIncrements": {   
      "OutputWriter" : [
        #{"format": "Paraview", "outputInterval": 1, "filename": "out/"+scenario_name+"/load_increments", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      ]
    },
  },
}
