import numpy as np

n = 4.
m = n

# boundary conditions (for quadratic elements)
bc = {}
for i in range(int(2*n+1)):
  x = i/(2*n+1.)
  bc[i] = np.sin(x*np.pi)
  i2 = (2*n+1)*(2*m) + i
  bc[i2] = np.sin(x*np.pi)
  
  print("bc[{}] = {}, bc[{}] = {}".format(i, bc[i], i2, bc[i2]))


# create megamol config file
config_file_contents = \
"""print('Hi, I am the megamolconfig.lua!')

-- mmSetAppDir("{megamol_home}/bin")
mmSetAppDir(".")

mmSetLogFile("")
mmSetLogLevel(0)
mmSetEchoLevel('*')

mmAddShaderDir("{megamol_home}/share/shaders")
mmAddResourceDir("{megamol_home}/share/resources")

mmPluginLoaderInfo("{megamol_home}/bin", "*.mmplg", "include")

-- mmSetConfigValue("*-window", "w1280h720")
mmSetConfigValue("*-window", "w720h720")
mmSetConfigValue("consolegui", "on")

mmSetConfigValue("LRHostEnable", "true")

return "done with megamolconfig.lua."
-- error("megamolconfig.lua is not happy!")
""".format(megamol_home="/store/software/opendihu/dependencies/megamol/install/install")

config_filename = "megamol_config.lua"
with open(config_filename, "w") as f:
  f.write(config_file_contents)

config = {
  "MegaMolArguments": "--configfile {} -p /store/software/opendihu/dependencies/megamol/src/megamol-master/project_files/testspheres.lua ".format(config_filename),
  "FiniteElementMethod" : {
    "nElements": [n, m],
    "inputMeshIsGlobal": True,
    "physicalExtent": [1.0, 1.0],
    "outputInterval": 1.0,
    "prefactor": 1,
    "dirichletBoundaryConditions": bc,
    "relativeTolerance": 1e-15,
    "solverType": "gmres",
    "preconditionerType": "none",
    "maxIterations": 10000,
    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/p", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":False},
      {"format": "PythonFile", "filename": "out/p", "outputInterval": 1, "binary":False, "onlyNodalValues":True},
      {"format": "MegaMol", "outputInterval": 1}
    ]
  },
}
