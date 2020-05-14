# isotropic Mooney Rivlin, using the FEBIO solver
import numpy as np
import sys, os

# parameters
force = 10
material_parameters = [10, 10, 1e5]       # c0, c1, k

# number of elements
nx = 2    # 2
ny = 2    # 2
nz = 5    # 5

physical_extent = [nx, ny, nz]

config = {
  "scenarioName":                   "3d_box_febio",
  "solverStructureDiagramFile":     "solver_structure.txt",             # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile":   "out/mappings_between_meshes.txt",  # output file for log of mappings 
  "NonlinearElasticitySolverFebio": {
    "durationLogKey": "febio",
    "force": force,                # factor of force that is applied in axial direction of the muscle
    "materialParameters": material_parameters,   # c0, c1, k for Ψ = c0 * (I1-3) + c1 * (I2-3) + 1/2*k*(log(J))^2
    
    # mesh
    "nElements": [nx, ny, nz],
    "inputMeshIsGlobal": True,
    "physicalExtent": physical_extent,
    "physicalOffset": [0, 0, 0],        # offset/translation where the whole mesh begins
    
    # output writer
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/febio", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      {"format": "PythonFile", "filename": "out/febio", "outputInterval": 1, "binary":False, "onlyNodalValues":True, "fileNumbering": "incremental"},
    ],
  },
}
