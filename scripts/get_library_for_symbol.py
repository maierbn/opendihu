#!/usr/bin/python3
#
# usage: get_library_for_symbol <link command> <symbol>
#
# Parses a compile or linker command and get the libraries from there. Find out in which library the given symbol is defined. The symbol can be either mangled or unmangled.
# An example link command is:
# g++ -o build_debug/solving_elasticity_problems_tutorial.cpp -pthread /store/software/opendihu/dependencies/easyloggingpp/install/src/easylogging++.o -fopenmp -Wl,-rpath=/store/software/opendihu/dependencies/lapack/install/lib -Wl,-rpath=/store/software/opendihu/dependencies/petsc/install/lib -Wl,-rpath=/store/software/opendihu/dependencies/python/install/lib -Wl,-rpath=/store/software/opendihu/dependencies/base64/install/lib -Wl,-rpath=/store/software/opendihu/dependencies/googletest/install/lib -Wl,-rpath=/store/software/opendihu/dependencies/semt/install/lib -Wl,-rpath=/store/software/opendihu/dependencies/adios/install/lib -Wl,-rpath=/store/software/opendihu/dependencies/vtk/install/lib -Wl,-rpath=/store/software/opendihu/dependencies/hdf5/install/lib -Wl,-rpath=/store/software/opendihu/dependencies/xercesc/install/lib -Wl,-rpath=/store/software/opendihu/dependencies/xsd/install/lib -Wl,-rpath=/store/software/opendihu/dependencies/boost/install/lib -Wl,-rpath=/store/software/opendihu/dependencies/chaste/install/lib/chaste -Wl,-rpath=/store/software/opendihu/dependencies/easyloggingpp/install/lib build_debug/src/solving_elasticity_problems_tutorial.o -L/usr/lib -L/usr/lib/x86_64-linux-gnu/openmpi/lib -L/store/software/opendihu/dependencies/lapack/install/lib -L/store/software/opendihu/dependencies/petsc/install/lib -L/store/software/opendihu/dependencies/python/install/lib -L/store/software/opendihu/dependencies/base64/install/lib -L/store/software/opendihu/dependencies/googletest/install/lib -L/store/software/opendihu/dependencies/semt/install/lib -L/store/software/opendihu/dependencies/adios/install/lib -L/store/software/opendihu/dependencies/vtk/install/lib -L/store/software/opendihu/dependencies/hdf5/install/lib -L/store/software/opendihu/dependencies/xercesc/install/lib -L/store/software/opendihu/dependencies/xsd/install/lib -L/store/software/opendihu/dependencies/boost/install/lib -L/store/software/opendihu/dependencies/chaste/install/lib/chaste -L/store/software/opendihu/dependencies/easyloggingpp/install/lib -L/store/software/opendihu/core/build_debug -lopendihud -lchaste_ode -lchaste_continuum_mechanics -lchaste_pde -lchaste_io -lchaste_mesh -lchaste_linalg -lchaste_global -lboost_wserialization -lboost_system -lboost_program_options -lboost_filesystem -lboost_serialization -lxerces-c -lvtkViewsInfovis-8.2 -lvtkChartsCore-8.2 -lvtkIOExportOpenGL2-8.2 -lvtkImagingStencil-8.2 -lvtkGeovisCore-8.2 -lvtkDomainsChemistryOpenGL2-8.2 -lvtkIOExportPDF-8.2 -lvtkIOExport-8.2 -lvtkRenderingContextOpenGL2-8.2 -lvtkViewsContext2D-8.2 -lvtkRenderingContext2D-8.2 -lvtkIOInfovis-8.2 -lvtkDomainsChemistry-8.2 -lvtkIOMINC-8.2 -lvtkIOAMR-8.2 -lvtkFiltersParallelImaging-8.2 -lvtkFiltersGeneric-8.2 -lvtkFiltersAMR-8.2 -lvtkIOParallel-8.2 -lvtkFiltersParallel-8.2 -lvtkIOAsynchronous-8.2 -lvtkInteractionImage-8.2 -lvtkIOCityGML-8.2 -lvtkIOExodus-8.2 -lvtkIOGeometry-8.2 -lvtkIOImport-8.2 -lvtkIOLSDyna-8.2 -lvtkIONetCDF-8.2 -lvtkIOParallelXML-8.2 -lvtkIOPLY-8.2 -lvtkIOSQL-8.2 -lvtkIOVideo-8.2 -lvtkRenderingGL2PSOpenGL2-8.2 -lvtkRenderingImage-8.2 -lvtkRenderingLabel-8.2 -lvtkRenderingLOD-8.2 -lvtkRenderingVolumeOpenGL2-8.2 -lvtkRenderingVolume-8.2 -lvtkFiltersPoints-8.2 -lvtkFiltersSMP-8.2 -lvtkFiltersTexture-8.2 -lvtkFiltersFlowPaths-8.2 -lvtkFiltersHyperTree-8.2 -lvtkFiltersImaging-8.2 -lvtkFiltersProgrammable-8.2 -lvtkFiltersSelection-8.2 -lvtkFiltersTopology-8.2 -lvtkFiltersVerdict-8.2 -lvtkImagingColor-8.2 -lvtkImagingMath-8.2 -lvtkImagingMorphological-8.2 -lvtkImagingStatistics-8.2 -lvtkIOEnSight-8.2 -lvtkIOMovie-8.2 -lvtkIOSegY-8.2 -lvtkIOTecplotTable-8.2 -lvtkIOVeraOut-8.2 -lvtkexodusII-8.2 -lvtkRenderingOpenGL2-8.2 -lvtkNetCDF-8.2 -lvtkParallelCore-8.2 -lvtkFiltersStatistics-8.2 -lvtkIOLegacy-8.2 -lvtkImagingFourier-8.2 -lvtkverdict-8.2 -lvtkgl2ps-8.2 -lvtkhdf5_hl-8.2 -lvtkhdf5-8.2 -lvtkIOXML-8.2 -lvtkViewsCore-8.2 -lvtkInteractionStyle-8.2 -lvtkInfovisLayout-8.2 -lvtkInteractionWidgets-8.2 -lvtkproj-8.2 -lvtkglew-8.2 -lvtkImagingSources-8.2 -lvtkImagingCore-8.2 -lvtkInfovisCore-8.2 -lvtkImagingHybrid-8.2 -lvtkFiltersExtraction-8.2 -lvtkRenderingAnnotation-8.2 -lvtkFiltersModeling-8.2 -lvtkFiltersHybrid-8.2 -lvtkImagingGeneral-8.2 -lvtkRenderingFreeType-8.2 -lvtkIOImage-8.2 -lvtkpugixml-8.2 -lvtklz4-8.2 -lvtklzma-8.2 -lvtkdoubleconversion-8.2 -lvtklibharu-8.2 -lvtkDICOMParser-8.2 -lvtkmetaio-8.2 -lvtkpng-8.2 -lvtktiff-8.2 -lvtkjpeg-8.2 -lvtklibxml2-8.2 -lvtktheora-8.2 -lvtkogg-8.2 -lvtkjsoncpp-8.2 -lvtksqlite-8.2 -lvtkIOXMLParser-8.2 -lvtkIOCore-8.2 -lvtkexpat-8.2 -lvtkzlib-8.2 -lvtkRenderingCore-8.2 -lvtkFiltersSources-8.2 -lvtkFiltersCore-8.2 -lvtkFiltersGeneral-8.2 -lvtkCommonColor-8.2 -lvtkFiltersGeometry-8.2 -lvtkCommonComputationalGeometry-8.2 -lvtkCommonCore-8.2 -lvtkCommonExecutionModel-8.2 -lvtkCommonSystem-8.2 -lvtkCommonDataModel-8.2 -lvtkCommonMath-8.2 -lvtkCommonTransforms-8.2 -lvtksys-8.2 -lvtkCommonMisc-8.2 -lvtkfreetype-8.2 -ladios2 -lgtest -lpython3.6m -lHYPRE -lzmumps -lmumps_common -lcmumps -lscalapack -lpord -lparmetis -ldmumps -lsundials_nvecparallel -lsundials_nvecserial -lpetsc -lsundials_cvode -lsmumps -lopenblas -lmpi_cxx -lmpi -lhdf5 -lpthread -ldl -lutil -lm


import sys, os
import subprocess

debug_output = False
verbose = False

if len(sys.argv) < 3:
  print("usage: {}  <link command> <symbol>".format(sys.argv[0]))
  sys.exit(0)

command = sys.argv[1]
symbol_to_search = bytes((sys.argv[2]).strip(), 'utf-8')

print("Symbol: \"{}\"".format(symbol_to_search))
if verbose:
  print("Command:")

libpaths = []
libnames = []

if os.environ.get("LD_LIBRARY_PATH") is not None:
  libpaths = os.environ.get("LD_LIBRARY_PATH").split(":")
  
libpaths.append("/usr/lib")

# parse command
for argument in command.split():
  argument = argument.strip()
  
  if verbose:
    print(argument)
  
  if "-L" == argument[0:2]:
    libpath = argument[2:]
    libpaths.append(libpath)
  elif "-l" == argument[0:2]:
    libname = argument[2:]
    libnames.append(libname)

if verbose:
  print("libpaths: {}".format(libpaths))
  print("libnames: {}".format(libnames))

  print("")
  print("Searching libraries")

# find libraries
libs = []
for i,libname in enumerate(libnames):
  lib = {
    "name": libname,
    "undefined_symbols": [],
    "known_symbols": [],
    "depending_on_libs": [],
    "depending_libs": [],
    "known_symbols_flags": [],
  }
  full_path = None
  for path in libpaths:
    full_path = os.path.join(path, "lib"+libname+".a")
    if os.path.isfile(full_path):
      lib["path"] = full_path
      lib["suffix"] = ".a"
      break
    full_path = os.path.join(path, "lib"+libname+".so")
    if os.path.isfile(full_path):
      lib["path"] = full_path
      lib["suffix"] = ".so"
      break
    full_path = None
  
  if full_path is None:
    if verbose:
      print("{}/{}: Could not find {}".format(i, len(libnames), libname))
    continue
  
  # get symbols
  try:
    output_mangled = subprocess.check_output("nm "+full_path, shell=True)
    output_unmangled = subprocess.check_output("nm -C "+full_path, shell=True)
    
    for output in [output_mangled, output_unmangled]:
      for line in output.split(bytes("\n", 'utf-8')):
        if ":" in str(line) or len(line) < 20:
          continue
        
        flag = str(line)[19]
        symbol = line[19:]
        
        
        #print("flag: {}, symbol: {}".format(flag, symbol))
        
        if flag == "U":
          lib["undefined_symbols"].append(symbol)
        else:
          lib["known_symbols"].append(symbol)
          lib["known_symbols_flags"].append(flag)
         
    if verbose:
      print("{}/{}: lib {} at path {}: {} known symbols, {} missing".format(i, len(libnames), libname, full_path, len(lib["undefined_symbols"]), len(lib["known_symbols"])))
    
    libs.append(lib)
  except:
    pass

print("Search for symbol \"{}\" in {} libraries.".format(symbol_to_search, len(libs)))

# find symbol
for lib in libs:
  
  if symbol_to_search in lib["known_symbols"]:
    for i,symbol in enumerate(lib["known_symbols"]):
      if symbol == symbol_to_search:
        break
    print("defined in {} at path {}, flag \"{}\"".format(lib["name"], lib["path"], lib["known_symbols_flags"][i]))
  
  if symbol_to_search in lib["undefined_symbols"]:
    print("needed by {} at path {}, flag \"U\"".format(lib["name"], lib["path"]))
