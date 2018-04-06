# 3D Incompressible Material
#
# command arguments: <name> <number elements>

import numpy as np
import scipy.integrate
import sys

n = 2   # number of elements
name = ""

if len(sys.argv) > 1:
  if "check_results.py" not in sys.argv[0]:
    name = sys.argv[0]
    n = int(sys.argv[1])
    
    print "name: \"{}\", n: {}".format(name, n)

nx = n
ny = n+1
nz = n+2

# for n=2: [2,3,4] = 24 elements, 8 dofs per element, 8*24=192 dofs, 192*3=576 unknows

nx = 1
ny = 1
nz = 1
# [1,1,1] = 1 element, 8 dofs per element, 8 dofs, 24 unknowns

config = {
  "ExplicitEuler" : {
    "initialValues": initial_values,
    "numberTimeSteps": 1500,
    "endTime": end_time,
    "FiniteElementMethod" : {
      "nElements": [nx,ny,nz],
      "nodeDimension": 1,
      "physicalExtent": [2.*nx,2.*ny,2.*nz],
      "DirichletBoundaryCondition": {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0},  # displacement Dirichlet bc
      "relativeTolerance": 1e-15,
      "rightHandSide": {0: -1.0},  # surface traction or body force
    },
    "OutputWriter" : [
      #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binaryOutput": "false", "fixedFormat": False},
      #{"format": "ExFile", "filename": "out/"+name, "outputInterval": 2},
      {"format": "PythonFile", "filename": "out/"+name, "outputInterval": 5, "binary":False, "onlyNodalValues":True},
    ]
  },
}

# output config in a readable format
if False:
  import pprint 
  pp = pprint.PrettyPrinter()
  pp.pprint(config)
