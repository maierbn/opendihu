
# isotropic Mooney Rivlin
import numpy as np
import sys, os

# parameters
force = 3
material_parameters = [2, 4]       # c0, c1

# number of elements
nx = 2    # 2
ny = 2    # 2
nz = 5    # 5

#physical_extent = [1, 1, 1]
physical_extent = [4, 4, 5]

# number of nodes
mx = 2*nx + 1
my = 2*ny + 1
mz = 2*nz + 1

# set the same Dirichlet boundary conditions as the FEBio settings
# boundary conditions (for quadratic elements)
dirichlet_bc = {}

xpos = 0.0
ypos = 0.0
zpos = 0.0

# fix z direction for all
for j in range(0,my):
  for i in range(0,mx):
    k = 0
    dirichlet_bc[k*mx*my + j*(mx) + i] = [None,None,0]

# fix x direction for left row
for j in range(0,my):
  k = 0
  i = 0
  dirichlet_bc[k*mx*my + j*mx + i][0] = 0

# fix y direction for front row
for i in range(0,mx):
  k = 0
  j = 0
  dirichlet_bc[k*mx*my + j*mx + i][1] = 0

# set the Neumann bc's
neumann_bc = [{"element": (nz-1)*nx*ny + j*nx + i, "constantVector": [0,0,force], "face": "2+"} for j in range(ny) for i in range(nx)]

def compare_result(data):
  
  import sys, os
  import py_reader
  import numpy as np

  #print(os.getcwd())
  # load data
  data1 = py_reader.load_data(["out/febio_0000001.py"])
  data2 = py_reader.load_data(["out/opendihu_0000001.py"])
  
  # if files do not exist
  if data1 == [] or data2 == []:
    return

  component_name = "0"
  total_error = 0

  values1 = py_reader.get_values(data1[0], "geometry", component_name)
  values2 = py_reader.get_values(data2[0], "geometry", component_name)
    
  # values2 contains entries for quadratic elements
  # extract the corner values
  n_elements = data2[0]['nElements']
  nx = n_elements[0]
  ny = n_elements[1]
  nz = n_elements[2]
  mx = nx*2 + 1
  my = ny*2 + 1
  mz = nz*2 + 1

  values2_linear = []
  for k in range(nz+1):
    for j in range(ny+1):
      for i in range(nx+1):
        values2_linear.append(values2[2*k * mx * my + 2*j * mx + 2*i])
    
  #print("values1 (febio):    ",list(values1))
  #print("values2 (opendihu): ",values2_linear)
    
  error_rms = np.sqrt(np.mean((values1-values2_linear)**2))
    
  print("rms: {}".format(error_rms))
  with open("rms", "w") as f:
    f.write(str(error_rms))

# ---------------------
# config for opendihu
config = {
  "scenarioName": "3d_box",
  "logFormat":                      "csv",     # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile":   "out/mappings_between_meshes.txt",  # output file for log of mappings 
  "HyperelasticitySolver": {
    "durationLogKey": "nonlinear",
    
    "materialParameters":         material_parameters,
    "displacementsScalingFactor": 1.0,   # scaling factor for displacements
    "constantBodyForce":          [0.0, 0.0, 0.0],   # body force in whole body region
    "residualNormLogFilename": "log_residual_norm.txt",
    "useAnalyticJacobian": True,
    "useNumericJacobian": False,   # Only works in parallel execution. If both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian
      
    "dumpDenseMatlabVariables": False,   # extra output of matlab vectors, x,r, jacobian matrix
    # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables are all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian
    
    # mesh
    "nElements": [nx, ny, nz],
    "inputMeshIsGlobal": True,
    "physicalExtent": physical_extent,
    "physicalOffset": [0, 0, 0],        # offset/translation where the whole mesh begins
    #"nodePositions": [[(i/(2.*nx))*(1+k/(2.*nz)), (j/(2*nx))*(1+k/(2.*nz)), (k/(2.*nz))**1.1*5] for k in range(mz) for j in range(my) for i in range(mx)],
    
    # nonlinear solver
    "relativeTolerance": 1e-10,         # 1e-10 relative tolerance of the linear solver
    "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual of the linear solver    
    "solverType": "preonly",            # type of the linear solver: cg groppcg pipecg pipecgrr cgne nash stcg gltr richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs fbcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly qcg bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres tsirm cgls
    "preconditionerType": "lu",         # type of the preconditioner
    "maxIterations": 1e4,               # maximum number of iterations in the linear solver
    "dumpFilename": "",
    "dumpFormat": "matlab",   # default, ascii, matlab
    "snesMaxFunctionEvaluations": 1e8,  # maximum number of function iterations
    "snesMaxIterations": 20,            # maximum number of iterations in the nonlinear solver
    "snesRelativeTolerance": 1e-5,     # relative tolerance of the nonlinear solver
    "snesLineSearchType": "l2",        # type of linesearch, possible values: "bt" "nleqerr" "basic" "l2" "cp" "ncglinear"
    "snesAbsoluteTolerance": 1e-5,     # absolute tolerance of the nonlinear solver
    "snesRebuildJacobianFrequency": 5, # how often the jacobian should be recomputed, -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time the Jacobian is built etc. -2 means rebuild at next chance but then never again 
    
    #"loadFactors":  [0.1, 0.2, 0.35, 0.47, 0.5, 0.75, 1.0],   # load factors for every timestep
    "loadFactors": [],                 # no load factors, solve problem directly
    "loadFactorGiveUpThreshold":    0.1,                # if the adaptive time stepping produces a load factor smaller than this value, the solution will be accepted for the current timestep, even if it did not converge fully to the tolerance
    "nNonlinearSolveCalls": 1,         # how often the nonlinear solve should be repeated
    
    # boundary conditions
    "dirichletBoundaryConditions": dirichlet_bc,
    "neumannBoundaryConditions": neumann_bc,
    "divideNeumannBoundaryConditionValuesByTotalArea": False,            # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
    #"updateDirichletBoundaryConditionsFunction": update_dirichlet_boundary_conditions,
    "updateDirichletBoundaryConditionsFunction": None,
    "updateDirichletBoundaryConditionsFunctionCallInterval": 1,
    
    "dirichletOutputFilename":     "out/dirichlet_boundary_conditions_tendon",    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
      
    "OutputWriter" : [   # output files for displacements function space (quadratic elements)
      {"format": "Paraview", "outputInterval": 1, "filename": "out/opendihu", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      {"format": "PythonFile", "filename": "out/opendihu", "outputInterval": 1, "binary":False, "onlyNodalValues":True, "fileNumbering": "incremental"},
      {"format": "PythonCallback", "callback": compare_result, "outputInterval": 1, "filename": "", "onlyNodalValues":True, "fileNumbering": "incremental"},
    ],
    "pressure": None,
    #"LoadIncrements": None,
    #"pressure": {   # output files for pressure function space (linear elements)
    #  "OutputWriter" : [
    #    {"format": "Paraview", "outputInterval": 1, "filename": "out/p", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
    #    {"format": "PythonFile", "filename": "out/p", "outputInterval": 1, "binary":False, "onlyNodalValues":True, "fileNumbering": "incremental"},
    #  ]
    #},
    # output writer for debugging, outputs files after each load increment, the geometry is not changed but u and v are written
    "LoadIncrements": {   
      "OutputWriter" : [
        {"format": "Paraview", "outputInterval": 1, "filename": "out/load_increments", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      ]
    },
  },
}

# ---------------------
# config for febio
import numpy as np
import sys, os

# parameters
force = 3
material_parameters = [2, 4, 1e6]       # c0, c1, k

# number of elements
nx = 2    # 2
ny = 2    # 2
nz = 5    # 5


config2 = {
  "NonlinearElasticitySolverFebio": {
    "durationLogKey": "febio",
    "force": force,                                     # factor of force that is applied in axial direction of the muscle
    "materialParameters": material_parameters,          # c0, c1, k for Ψ = c0 * (I1-3) + c1 * (I2-3) + 1/2*k*(log(J))^2
    "tractionElementNos": [(nz-1)*nx*ny + j*nx + i for j in range(ny) for i in range(nx)],    # elements on which traction is applied
    "tractionVector": [0,0,force],                      # traction vector that is applied
    "dirichletBoundaryConditionsMode": "fix_floating",  # "fix_all" or "fix_floating", how the bottom of the box will be fixed, fix_all fixes all nodes, fix_floating fixes all nodes only in z and the edges in x/y direction
    "slotNames": [],
    
    # mesh
    "nElements": [nx, ny, nz],
    "inputMeshIsGlobal": True,
    "physicalExtent": physical_extent,
    "physicalOffset": [0, 0, 0],        # offset/translation where the whole mesh begins
    #"nodePositions": [[i/nx*(1+k/nz), j/ny*(1+k/nz), (k/nz)**1.1*5] for k in range(nz+1) for j in range(ny+1) for i in range(nx+1)],
    
    # output writer
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/febio", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      {"format": "PythonFile", "filename": "out/febio", "outputInterval": 1, "binary":False, "onlyNodalValues":True, "fileNumbering": "incremental"},
    ],
  },
}

config.update(config2)
