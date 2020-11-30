#!/bin/bash
#
# This scripts creates fibers as streamlines in a potential flow simulation.
# Input is an STL mesh of the surface of the muscle. Output is the pickle file that contains the fibers.
# This is not the algorithm described in the dissertation of BM and the result is not a 3D mesh, but individual 1D fiber meshes with different lengths.
#
# This file needs to be run from the directory "examples/fiber_tracing/streamline_tracer/scripts"

# directories
basedir=$(pwd)/..
opendihu_directory=$(pwd)/../../../..
stl_utility_directory=${opendihu_directory}/scripts/stl_utility
pyod=${opendihu_directory}/dependencies/python/install/bin/python3
# -----------------
# settings

# z range of the mesh for which fibers should be extracted
min_z=40
max_z=260

# number of subdivisions of the z range, this has to be odd for quadratic elements
n_rings=13
#n_rings=101
#n_rings=43

# choose which algorithm to use
triangulation_type=2         # triangulation_type:  0 = scipy, 1 = triangle, 2 = custom with CoG, 3 = custom with minimized distance
parametric_space_shape=3     # 0 = unit circle, 1 = unit square, 2 = unit square with adjusted grid, 3 = unit circle with adjusted grid, 4 = unit circle with optimized positions

# number of fibers in x or y direction,  has to be even for quadratic elements
#n_points_x=100
#n_points_x=32
#n_points_x=10
n_points_x=4

# input mesh to use, this is an STL mesh of the surface of the geometry
input_file=${basedir}/original_meshes/biceps_full.stl

# prefix for output file
prefix="biceps_"

# use the following to create it:
#$pyod ${stl_utility_directory}/remove_inside_triangles.py ${input_file} ${input_file}_repaired.stl    # remove non-surface triangles (takes ~1h)
# -----------------

# change to working directory for meshes
mkdir -p ${basedir}/processed_meshes/out
mkdir -p ${basedir}/processed_meshes/results
cd ${basedir}/processed_meshes

cp $input_file mesh_00.stl    # cp input file

# ---- python part, preprocessing ------
# cut surface mesh at specified z positions and create rings from it. Write result to `rings_created`, debugging output to out/mesh_01.stl (takes <1min)
# usage: ./create_rings.py <input file> [<n rings> [<min z> <max z>]]
$pyod ${basedir}/scripts/utility/create_rings.py $input_file $n_rings $min_z $max_z

# extract the existing rings from the surface mesh, i.e. use the nodes in the STL mesh
# (this is the alternative to create_rings.py)   
# write result to `rings_extracted`.
#$pyod ${basedir}/scripts/utility/extract_rings.py ../biceps.stl     

# rename ring output file to `rings`
mv rings_created rings      
#mv rings_extracted rings

# create a mesh, reads input from `rings`, write output to `mesh`, debugging output to out/mesh_02* - out/mesh_09* (takes ~2min)
# usage: ./create_mesh.py [<triangulation_type> [<parametric_space_shape> [<n_points_x> [<n_grid_points_x>]]]]"
$pyod ${basedir}/scripts/utility/create_mesh.py $triangulation_type $parametric_space_shape $n_points_x

input_mesh_name="${basedir}/processed_meshes/mesh"

# ---- c++ part, streamline tracing ------
#ln -s mesh mesh_normal

# compile program
cd $opendihu_directory
scons no_tests=TRUE

# change to basedir and compile example
cd ${basedir}
scons
cd build_release

# run the simulation with different discretization orders and structured or unstructured meshes, they produce all more or less the same output
# the following is a list of the executables
schemes="\
laplace_structured_linear \
laplace_structured_quadratic \
laplace_structured_hermite \
laplace_unstructured_linear \
laplace_unstructured_quadratic \
laplace_unstructured_hermite"

for scheme_name in $schemes; do

  # run the simulation, created meshes with 100. elements per cm and the longest with 15 cm length

  # command arguments: <input_filename> <output_filename> [<target_element_length> [<target_fiber_length>]]
  echo ./${scheme_name} ../settings_streamline_tracer.py ${input_mesh_name} ${prefix}${scheme_name} 1e-2 15
  ./${scheme_name} ../settings_streamline_tracer.py ${input_mesh_name} ${prefix}${scheme_name} 1e-2 15

  # convert csv output file to stl and pickle file
  $pyod ${basedir}/scripts/utility/convert_mesh_csv_to_stl.py ${prefix}${scheme_name}.csv ../processed_meshes/results/${prefix}${scheme_name}.stl
  $pyod ${basedir}/scripts/utility/convert_mesh_csv_to_pickle.py ${prefix}${scheme_name}.csv ../processed_meshes/results/${prefix}${scheme_name}.pickle
  $pyod ${basedir}/scripts/utility/convert_mesh_csv_to_pickle.py ${prefix}${scheme_name}_raw.csv ../processed_meshes/results/${prefix}${scheme_name}_raw.pickle

done



