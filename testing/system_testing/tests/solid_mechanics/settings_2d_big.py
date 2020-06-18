# 2D Incompressible Material
#
# command arguments: <analytical jacobian> <name>

import numpy as np
import scipy.integrate
import sys,os

name = ""
analytic_jacobian = True

if len(sys.argv) > 1:
  if "check_results.py" in sys.argv[0]:
    del sys.argv[0]
      
if len(sys.argv) > 1:
  analytic_jacobian = True if int(sys.argv[0]) == 1 else False
  name = os.path.basename(sys.argv[1])
  
  print("analytic jacobian: {}".format(analytic_jacobian))
  print("name: \"{}\"".format(name))
    
nx = 10
ny = 4
# [10,4] = 40 elements, 55 nodes, 1 dof per node, 2 unknowns per dof, 110 unknowns

lx = 8.0
ly = 3.0
tmax = 0.5
that = 0.5

dirichletBC = {
  0: 0.0, 1: 0.0,
  22: 0.0,
  44: 0.0,
  66: 0.0,
  88: 0.0,
} 

material_parameters = [1.0, 0.0, 1.0]  # c0, c1, kappa
# kappa = 0 means disable incompressibility and volume effects

config = {
  "logFormat": "csv",
  "FiniteElementMethod" : {
    "nElements": [nx,ny],
    "physicalExtent": [lx,ly],
    "dirichletBoundaryCondition": dirichletBC,  # displacement Dirichlet bc
    "tractionReferenceConfiguration": [
      {"element": 9, "face": "0+", "constantValue": tmax/(ly/4.)},   # face can be one of "0+", "0-", "1+", "1-", "2+", "2-". dofValues are the element-local node/dof numbers (e.g. 0-26 for 3D quadratic Lagrange elements), the entries are vectors of dimension D with respect to the local element coordinate system
      {"element": 19, "face": "0+", "constantValue": tmax/(ly/4.)},   # face can be one of "0+", "0-", "1+", "1-", "2+", "2-". dofValues are the element-local node/dof numbers (e.g. 0-26 for 3D quadratic Lagrange elements), the entries are vectors of dimension D with respect to the local element coordinate system
      {"element": 29, "face": "0+", "constantValue": tmax/(ly/4.)},   # face can be one of "0+", "0-", "1+", "1-", "2+", "2-". dofValues are the element-local node/dof numbers (e.g. 0-26 for 3D quadratic Lagrange elements), the entries are vectors of dimension D with respect to the local element coordinate system
      {"element": 39, "face": "0+", "constantValue": tmax/(ly/4.)},   # face can be one of "0+", "0-", "1+", "1-", "2+", "2-". dofValues are the element-local node/dof numbers (e.g. 0-26 for 3D quadratic Lagrange elements), the entries are vectors of dimension D with respect to the local element coordinate system
    ],
    #"tractionCurrentConfiguration": [
    #  {"element": 0, "face": "0+", "constantValue": that}   # face can be one of "0+", "0-", "1+", "1-", "2+", "2-". dofValues are the element-local node/dof numbers (e.g. 0-26 for 3D quadratic Lagrange elements), the entries are vectors of dimension D with respect to the local element coordinate system
    #],
    #"tractionReferenceConfiguration": [
    #  {"element": 1, "face": "0+", "dofVectors": {0: [tmax,0,0], 1: [tmax,0,0], 2: [tmax,0,0], 3: [tmax,0,0]}},  # face can be one of "0+", "0-", "1+", "1-", "2+", "2-". dofValues are the element-local node/dof numbers (e.g. 0-26 for 3D quadratic Lagrange elements), the entries are vectors of dimension D with respect to the local element coordinate system
    #  {"element": 1, "face": "1-", "dofVectors": {5: [1,2,3], 7:[8,3,4]}}, 
    #  {"element": 3, "face": "0-", "constantVector": tmax},   # a constant vector for constant traction on that element
    #  {"element": 3, "face": "0+", "constantValue": 5.12},   # a single constant value for constant traction on that element, in direction of outward normal
    #],
    #"tractionCurrentConfiguration": [
    #  
    #],
    #"bodyForceReferenceConfiguration": {0: [tmax,0,0], 5: [tmax,tmax,tmax]},   # {<element global no.>: <force vector, constant in element>, ...}
    #"bodyForceCurrentConfiguration": {},
    "relativeTolerance": 1e-15,
    #"rightHandSide": traction,  # surface traction T or body force B, both in material description
    "materialParameters": material_parameters,  # c0, c1, kappa
    "analyticJacobian": analytic_jacobian,   # False = compute Jacobian by finite differences
      
    "OutputWriter" : [
      #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binary": "false", "fixedFormat": False},
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
