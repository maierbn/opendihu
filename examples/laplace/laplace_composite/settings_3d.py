import numpy as np

# boundary conditions (for quadratic elements)
bc = {}
n_nodes_x1 = 5
for i in range(int(n_nodes_x1)):
  x = i/n_nodes_x1
  bc[i] = np.sin(x*np.pi)*-3
  
# z+ top
for j in range(3):
  for i in range(5):
    bc[15*6+j*5+i] = 2.0
  
# y+ top
for i in range(5):
  bc[15*7+7+i] = -1.0
  
# x+ top
for k in range(5):
  bc[15*7+k*7*3 + 7*3-1] = 1.0
  bc[15*7+k*7*3 + 7*3-2] = 1.0
  
  
config = {
  "solverStructureDiagramFile":     "solver_structure.txt",     # output file of a diagram that shows data connection between solvers
  "logFormat":                      "csv",                      # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "scenarioName":                   "laplace_composite",        # scenario name to find the run in the log file
  "mappingsBetweenMeshesLogFile":   "",                         # a log file about mappings between meshes, here we do not want that because there are no mappings
  "Meshes": {
    "submesh0": {
      "nElements": [2,1,3],
      "inputMeshIsGlobal": True,
      "physicalExtent": [4,2,6],
      "physicalOffset": [0,0,0],
    },
    "submesh1": {
      "nElements": [3,1,2],
      "inputMeshIsGlobal": True,
      "nodePositions": 
        sum([[[x,2,z] for x in range(5)] + [[4,1,z],[4,0,z]]
        + [[x,3,z] for x in range(4)] + [[5,2.5,z],[5,1,z],[5,0,z]]
        + [[x,4,z] for x in range(4)] + [[6,3,z],[6,1,z],[6,0,z]] for z in range(5)], [])
    },
  },
  "FiniteElementMethod": {
    "inputMeshIsGlobal": True,
    "outputInterval": 1.0,
    
    "dirichletBoundaryConditions": bc,
    "dirichletOutputFilename":     None,                                # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
    "neumannBoundaryConditions": [],
    "prefactor": [1,2],
    
    "meshName": ["submesh0", "submesh1"],
    
    # solver parameters
    "solverType": "gmres",
    "preconditionerType": "none",
    "relativeTolerance": 1e-15,
    "absoluteTolerance": 1e-10,         # 1e-10 absolute tolerance of the residual        
    "maxIterations": 10000,
    "dumpFormat": "default",
    "dumpFilename": "",
    "slotName": "",
    
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/3d", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      #{"format": "PythonFile", "filename": "out/laplace", "outputInterval": 1, "binary":False, "onlyNodalValues":True, "fileNumbering": "incremental"},
    ]
  },
}
