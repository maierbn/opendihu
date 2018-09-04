# Diffusion 2D Debug
n = 20   # number of elements

# initial values
iv = {}

for y in range(int(0.2*n), int(0.3*n)):
  for x in range(int(0.5*n), int(0.8*n)):
    i = y*(n+1) + x
    iv[i] = 1.0

print("iv: ",iv)

config = {
  "Heun" : {
    "initialValues": iv,
    "timeStepWidth": 1e-3,
    "endTime": 1.0,
    
    "FiniteElementMethod" : {
      "nElements": [n,n],
      "physicalExtent": [4.0,4.0],
      "relativeTolerance": 1e-15,
      "prefactor": 0.1,
    },
    
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/paraview/diffusion2d", "binary": "false", "fixedFormat": False},
      {"format": "PythonFile", "outputInterval": 10, "filename": "out/python/out_diffusion2d", "binary": True}
    ]
  },
}
