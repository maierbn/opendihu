nx = 10
ny = 5
nz = 3

dirichletBC = {
  0: 0.0, 1: 0.0,
  6: 0.0,
  12: 0.0,
  18: 0.0,
  24: 0.0,
  30: 0.0,
  36: 0.0,
  42: 0.0,
  48: 0.0,
} 

material_parameters = [5.0, 3.0, 100.0]  # c0, c1, kappa

config = {
  "logFormat": "csv",
  "FiniteElementMethod" : {
    "nElements": [nx,ny,nz],
    "physicalExtent": [10.0,10.0,10.0],
    "dirichletBoundaryCondition": dirichletBC,  # displacement Dirichlet bc
    "tractionReferenceConfiguration": {"element": 0, "face": "0+", "constantValue": 0.1},
    "relativeTolerance": 1e-15,
    "materialParameters": material_parameters,  # c0, c1, kappa
    "analyticJacobian": True,
    "numericJacobian": True,
    "logfile": "residual_norm.txt",
      
    "OutputWriter" : [
      {"format": "PythonFile", "filename": "out/nonlinear_test", "outputInterval": 1, "binary":False, "onlyNodalValues":True},
    ],
    "outputIntermediateSteps": True
  },
}
