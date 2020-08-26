# Linear elasticity
# Use paraview and the Warp Filter for visualization.
# Add a Glyph filter with arrows for field "-rhsNeumannBC", this is the applied traction.

import numpy as np

nx = 4
ny = 4

# boundary conditions (for linear elements)
dirichlet_bc = {0: [0.0,0.0,0.0]}

for j in range(1,ny+1):
  dirichlet_bc[j*(nx+1)] = [0.0,None,None]

dirichlet_bc[1] = [0.0,0.0]

neumann_bc = [{"element": j*nx+(nx-1), "constantVector": [0.1,+0.2], "face": "0+"} for j in range(ny)]

#dirichlet_bc = {}
#neumann_bc = []

config = {
  "scenarioName":                   "linear_elasticity_2d",    # scenario name for the log file
  "logFormat":                      "csv",     # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile":   "mappings_between_meshes.txt",   # log file for mappings between meshes
  "FiniteElementMethod": {
    
    # mesh
    "nElements":          [nx, ny],                 # number of elements in the 2 coordinate directions  
    "inputMeshIsGlobal":  True,                     # if the specifiation of the mesh is given in global numbers
    "physicalExtent":     [nx, ny],                 # the physical dimensions of the mesh, in this case the same as the number of elements
    "physicalOffset":     [0, 0, 0],                # lower left bottom coordinates of the mesh
    "outputInterval":     1.0,                      # the timestep interval in which to print the current time step to output
    "prefactor": 1,                                 # prefactor c of the equation c*Δu = f
    "dirichletBoundaryConditions":  dirichlet_bc,   # dirichlet boundary conditions
    "dirichletOutputFilename":      None,           # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    "neumannBoundaryConditions":    neumann_bc,     # neumann boundary conditions
    "divideNeumannBoundaryConditionValuesByTotalArea": False,  # if the neumann boundary condition vectors should be divided by the total surface area where surface loads are applied, this allows to specify the total force that acts on the surface. If set to False (default), the given traction is a per-surface quantity.
    "slotName":           "",                       # slot name of the solution variable
    
    # solver
    "relativeTolerance":  1e-15,                    # relative tolerance for the solver, this relates to the current residual norm divided by the norm of the rhs
    "absoluteTolerance":  1e-10,                    # absolute tolerance of the residual        
    "solverType":         "gmres",                  # which PETSc linear system solver to use
    "preconditionerType": "none",                   # which PETSc preconditioner to use
    "maxIterations":      1e4,                      # maximum number of iterations of the linear solver
    "dumpFormat":         "matlab",                 # "default", "ascii", "matlab", output format type for system matrix and right hand side
    "dumpFilename":       "",                       # filename for output of system matrix and right hand side after every solve, set to "" to disable
    
    # material parameters
    "bulkModulus":        1.5,                      # bulk modulus, K, material parameter for compressibility
    "shearModulus":       2.0,                      # shear modulus, μ
    
    "OutputWriter": [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/2d", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      {"format": "PythonFile", "filename": "out/2d", "outputInterval": 1, "binary":False, "onlyNodalValues":True, "fileNumbering": "incremental"},
    ]
  },
}
