#!/usr/bin/python3
#
# usage: get_static_link_order <link command>
#
# Parses a compile or linker command and rearranges the statically linked libraries such that dependencies are met.
# An example link command is:
# g++ -o build_debug/solving_elasticity_problems_tutorial.cpp -pthread /store/software/opendihu/dependencies/easyloggingpp/install/src/easylogging++.o -fopenmp -Wl,-rpath=/store/software/opendihu/dependencies/lapack/install/lib -Wl,-rpath=/store/software/opendihu/dependencies/petsc/install/lib -Wl,-rpath=/store/software/opendihu/dependencies/python/install/lib -Wl,-rpath=/store/software/opendihu/dependencies/base64/install/lib -Wl,-rpath=/store/software/opendihu/dependencies/googletest/install/lib -Wl,-rpath=/store/software/opendihu/dependencies/semt/install/lib -Wl,-rpath=/store/software/opendihu/dependencies/adios/install/lib -Wl,-rpath=/store/software/opendihu/dependencies/vtk/install/lib -Wl,-rpath=/store/software/opendihu/dependencies/hdf5/install/lib -Wl,-rpath=/store/software/opendihu/dependencies/xercesc/install/lib -Wl,-rpath=/store/software/opendihu/dependencies/xsd/install/lib -Wl,-rpath=/store/software/opendihu/dependencies/boost/install/lib -Wl,-rpath=/store/software/opendihu/dependencies/chaste/install/lib/chaste -Wl,-rpath=/store/software/opendihu/dependencies/easyloggingpp/install/lib build_debug/src/solving_elasticity_problems_tutorial.o -L/usr/lib -L/usr/lib/x86_64-linux-gnu/openmpi/lib -L/store/software/opendihu/dependencies/lapack/install/lib -L/store/software/opendihu/dependencies/petsc/install/lib -L/store/software/opendihu/dependencies/python/install/lib -L/store/software/opendihu/dependencies/base64/install/lib -L/store/software/opendihu/dependencies/googletest/install/lib -L/store/software/opendihu/dependencies/semt/install/lib -L/store/software/opendihu/dependencies/adios/install/lib -L/store/software/opendihu/dependencies/vtk/install/lib -L/store/software/opendihu/dependencies/hdf5/install/lib -L/store/software/opendihu/dependencies/xercesc/install/lib -L/store/software/opendihu/dependencies/xsd/install/lib -L/store/software/opendihu/dependencies/boost/install/lib -L/store/software/opendihu/dependencies/chaste/install/lib/chaste -L/store/software/opendihu/dependencies/easyloggingpp/install/lib -L/store/software/opendihu/core/build_debug -lopendihud -lchaste_ode -lchaste_continuum_mechanics -lchaste_pde -lchaste_io -lchaste_mesh -lchaste_linalg -lchaste_global -lboost_wserialization -lboost_system -lboost_program_options -lboost_filesystem -lboost_serialization -lxerces-c -lvtkViewsInfovis-8.2 -lvtkChartsCore-8.2 -lvtkIOExportOpenGL2-8.2 -lvtkImagingStencil-8.2 -lvtkGeovisCore-8.2 -lvtkDomainsChemistryOpenGL2-8.2 -lvtkIOExportPDF-8.2 -lvtkIOExport-8.2 -lvtkRenderingContextOpenGL2-8.2 -lvtkViewsContext2D-8.2 -lvtkRenderingContext2D-8.2 -lvtkIOInfovis-8.2 -lvtkDomainsChemistry-8.2 -lvtkIOMINC-8.2 -lvtkIOAMR-8.2 -lvtkFiltersParallelImaging-8.2 -lvtkFiltersGeneric-8.2 -lvtkFiltersAMR-8.2 -lvtkIOParallel-8.2 -lvtkFiltersParallel-8.2 -lvtkIOAsynchronous-8.2 -lvtkInteractionImage-8.2 -lvtkIOCityGML-8.2 -lvtkIOExodus-8.2 -lvtkIOGeometry-8.2 -lvtkIOImport-8.2 -lvtkIOLSDyna-8.2 -lvtkIONetCDF-8.2 -lvtkIOParallelXML-8.2 -lvtkIOPLY-8.2 -lvtkIOSQL-8.2 -lvtkIOVideo-8.2 -lvtkRenderingGL2PSOpenGL2-8.2 -lvtkRenderingImage-8.2 -lvtkRenderingLabel-8.2 -lvtkRenderingLOD-8.2 -lvtkRenderingVolumeOpenGL2-8.2 -lvtkRenderingVolume-8.2 -lvtkFiltersPoints-8.2 -lvtkFiltersSMP-8.2 -lvtkFiltersTexture-8.2 -lvtkFiltersFlowPaths-8.2 -lvtkFiltersHyperTree-8.2 -lvtkFiltersImaging-8.2 -lvtkFiltersProgrammable-8.2 -lvtkFiltersSelection-8.2 -lvtkFiltersTopology-8.2 -lvtkFiltersVerdict-8.2 -lvtkImagingColor-8.2 -lvtkImagingMath-8.2 -lvtkImagingMorphological-8.2 -lvtkImagingStatistics-8.2 -lvtkIOEnSight-8.2 -lvtkIOMovie-8.2 -lvtkIOSegY-8.2 -lvtkIOTecplotTable-8.2 -lvtkIOVeraOut-8.2 -lvtkexodusII-8.2 -lvtkRenderingOpenGL2-8.2 -lvtkNetCDF-8.2 -lvtkParallelCore-8.2 -lvtkFiltersStatistics-8.2 -lvtkIOLegacy-8.2 -lvtkImagingFourier-8.2 -lvtkverdict-8.2 -lvtkgl2ps-8.2 -lvtkhdf5_hl-8.2 -lvtkhdf5-8.2 -lvtkIOXML-8.2 -lvtkViewsCore-8.2 -lvtkInteractionStyle-8.2 -lvtkInfovisLayout-8.2 -lvtkInteractionWidgets-8.2 -lvtkproj-8.2 -lvtkglew-8.2 -lvtkImagingSources-8.2 -lvtkImagingCore-8.2 -lvtkInfovisCore-8.2 -lvtkImagingHybrid-8.2 -lvtkFiltersExtraction-8.2 -lvtkRenderingAnnotation-8.2 -lvtkFiltersModeling-8.2 -lvtkFiltersHybrid-8.2 -lvtkImagingGeneral-8.2 -lvtkRenderingFreeType-8.2 -lvtkIOImage-8.2 -lvtkpugixml-8.2 -lvtklz4-8.2 -lvtklzma-8.2 -lvtkdoubleconversion-8.2 -lvtklibharu-8.2 -lvtkDICOMParser-8.2 -lvtkmetaio-8.2 -lvtkpng-8.2 -lvtktiff-8.2 -lvtkjpeg-8.2 -lvtklibxml2-8.2 -lvtktheora-8.2 -lvtkogg-8.2 -lvtkjsoncpp-8.2 -lvtksqlite-8.2 -lvtkIOXMLParser-8.2 -lvtkIOCore-8.2 -lvtkexpat-8.2 -lvtkzlib-8.2 -lvtkRenderingCore-8.2 -lvtkFiltersSources-8.2 -lvtkFiltersCore-8.2 -lvtkFiltersGeneral-8.2 -lvtkCommonColor-8.2 -lvtkFiltersGeometry-8.2 -lvtkCommonComputationalGeometry-8.2 -lvtkCommonCore-8.2 -lvtkCommonExecutionModel-8.2 -lvtkCommonSystem-8.2 -lvtkCommonDataModel-8.2 -lvtkCommonMath-8.2 -lvtkCommonTransforms-8.2 -lvtksys-8.2 -lvtkCommonMisc-8.2 -lvtkfreetype-8.2 -ladios2 -lgtest -lpython3.6m -lHYPRE -lzmumps -lmumps_common -lcmumps -lscalapack -lpord -lparmetis -ldmumps -lsundials_nvecparallel -lsundials_nvecserial -lpetsc -lsundials_cvode -lsmumps -lopenblas -lmpi_cxx -lmpi -lhdf5 -lpthread -ldl -lutil -lm

import sys, os
import subprocess

debug_output = False

if len(sys.argv) < 2:
  print("usage: {} <link command>".format(sys.argv[0]))
  sys.exit(0)

command = sys.argv[1]

libpaths = []
libnames = []

if os.environ.get("LD_LIBRARY_PATH") is not None:
  libpaths = os.environ.get("LD_LIBRARY_PATH").split(":")
  
libpaths.append("/usr/lib")

# parse command
for argument in command.split():
  argument = argument.strip()
  
  print(argument)
  
  if "-L" == argument[0:2]:
    libpath = argument[2:]
    libpaths.append(libpath)
  elif "-l" == argument[0:2]:
    libname = argument[2:]
    libnames.append(libname)

print("libpaths: {}".format(libpaths))
print("libnames: {}".format(libnames))

# find libraries
libs = []
for i,libname in enumerate(libnames):
  lib = {
    "name": libname,
    "undefined_symbols": [],
    "known_symbols": [],
    "depending_on_libs": [],
    "depending_libs": [],
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
    print("{}/{}: Could not find {}".format(i, len(libnames), libname))
    continue
  
  # get symbols
  output = subprocess.check_output("nm "+full_path, shell=True)
  for line in output.split(bytes("\n", 'utf-8')):
    if ":" in str(line) or len(line) < 20:
      continue
    
    flag = str(line)[19]
    symbol = line[19:]
    
    #print("flag: {}, symbol: {}".format(flag, symbol))
    
    if flag == "U":
      lib["undefined_symbols"].append(symbol)
    elif flag == "T" or flag == "t":
      lib["known_symbols"].append(symbol)
   
  print("{}/{}: lib {} at path {}: {} known symbols, {} missing".format(i, len(libnames), libname, full_path, len(lib["undefined_symbols"]), len(lib["known_symbols"])))
  
  libs.append(lib)

# get dependend libs
for i,lib in enumerate(libs):
  for undefined_symbol in lib["undefined_symbols"]:
    for j,lib2 in enumerate(libs):
      if lib2["name"] == lib["name"] or lib2["name"] in libs[i]["depending_on_libs"]:
        continue
      if undefined_symbol in lib2["known_symbols"]:
        libs[i]["depending_on_libs"].append(lib2["name"])
        libs[j]["depending_libs"].append(lib["name"])
        print("  {} depends on {}".format(lib["name"], lib2["name"]))
        break
  
# print result
if debug_output:
  for lib in libs:
    print("{} before {}".format(lib["name"], lib["depending_on_libs"]))

# determine order
ordered = [[]]
for i,lib in enumerate(libs):
  
  if debug_output:
    print("{}/{} {} before {}".format(i, len(libs), lib["name"], lib["depending_on_libs"]))
  
  already_in_ordered = False
  for j,bucket in enumerate(ordered):
    if lib["name"] in bucket:
      already_in_ordered = True
      
      if debug_output:
        print("already in ordered in bucket {} {}".format(j, bucket))
        
      bucket_split_2nd_part = []
      to_append = []
      for depending_on_libs in lib["depending_on_libs"]:
        if depending_on_libs in bucket:
          bucket_split_2nd_part.append(depending_on_libs)
        else:
          to_append.append(depending_on_libs)
          
      bucket_split_1st_part = []
      for b in bucket:
        if b not in bucket_split_2nd_part:
          bucket_split_1st_part.append(b)
            
      append_at_end = []
      for l in to_append:
        l_is_already_in_later_bucket = False
        for bucket2 in ordered[j+1:]:
          if l in bucket2:
            l_is_already_in_later_bucket = True
            break
        if not l_is_already_in_later_bucket:
          append_at_end.append(l)
          
      if debug_output:
        print("add {} to bucket {} -> {},{}, other: {}, append: {}".format(lib["name"], bucket, bucket_split_1st_part, bucket_split_2nd_part, to_append, append_at_end))
        
      new_ordered = ordered[0:j] + [bucket_split_1st_part]
      if bucket_split_2nd_part:
        new_ordered += [bucket_split_2nd_part]
      new_ordered += ordered[j+1:]
      if append_at_end:
        new_ordered += [append_at_end]
      ordered = new_ordered
  
  if not already_in_ordered:
    if debug_output:
      print("not already in ordered")
    
    # find first bucket where a depending_on_libs is
    for j,bucket in enumerate(ordered):
      
      a_depending_lib_is_in_this_bucket = False
      for depending_on_libs in lib["depending_on_libs"]:
        if depending_on_libs in bucket:
          a_depending_lib_is_in_this_bucket = True
          
          if debug_output:
            print("a_depending_lib_is_in_this_bucket, lib {} in bucket {}".format(depending_on_libs, bucket))
          break
      if a_depending_lib_is_in_this_bucket:
        break
    
    if debug_output:
      print("ordered: {}, bucket {} of {}".format(ordered, j,len(ordered)))
    if len(ordered) == 0:
      j = 0
    elif j >= len(ordered):
      j = len(ordered)-1
      bucket = ordered[j]
      
    if debug_output:
      print("bucket: {}".format(bucket))
    
    # add to this bucket
    bucket_split_2nd_part = []
    to_append = []
    for depending_on_libs in lib["depending_on_libs"]:
      if depending_on_libs in bucket:
        bucket_split_2nd_part.append(depending_on_libs)
      else:
        to_append.append(depending_on_libs)
        
    
    bucket_split_1st_part = [lib["name"]]
    for b in bucket:
      if b not in bucket_split_2nd_part:
        bucket_split_1st_part.append(b)
        
    append_at_end = []
    for l in to_append:
      l_is_already_in_later_bucket = False
      for bucket2 in ordered[j+1:]:
        if l in bucket2:
          l_is_already_in_later_bucket = True
          break
      if not l_is_already_in_later_bucket:
        append_at_end.append(l)
        
    if debug_output:
      print("add {} to bucket {} -> {},{}, other: {}, append: {}".format(lib["name"], bucket, bucket_split_1st_part, bucket_split_2nd_part, to_append, append_at_end))
        
    new_ordered = ordered[0:j] + [bucket_split_1st_part]
    if bucket_split_2nd_part:
      new_ordered += [bucket_split_2nd_part]
    new_ordered += ordered[j+1:]
    if append_at_end:
      new_ordered += [append_at_end]
    ordered = new_ordered
    
  if debug_output:
    print("after {}: {}".format(lib["name"], ordered))
        
print("final order, items in each line have no defined order:")
final_order = []
for bucket in ordered:
  print(bucket)
  for lib in bucket:
    final_order.append(lib)
print("")
print("as list:")
print(final_order)
print("")

  
resulting_command = ""
libraries_appended = False
for argument in command.split(" "):
  argument = argument.strip()
  
  if "-l" == argument[0:2] and not libraries_appended:
    for bucket in ordered:
      for lib in bucket:
        resulting_command = resulting_command + "-l" + lib + " "
    libraries_appended = True
  else:
    resulting_command = resulting_command + str(argument) + " "
    
print("resulting command:")
print("")
print(resulting_command)

