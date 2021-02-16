#!/bin/bash
#
# This scripts creates fibers as streamlines in a potential flow simulation.
# Input is an STL mesh of the surface of the muscle. Output is the pickle file that contains the fibers.
#
# arguments <input_stl_file> <output_pickle_file> <output_bin_file> <min_z> <max_z> <n_points_x> <n_points_z> [--only-stage-1]
#
# if --only-stage-1 is set, only the 3D mesh will be created and no streamlines will be traced.

if [ "$#" -lt "6" ]; then
  echo "$# usage: <input_stl_file> <output_pickle_file> <min_z> <max_z> <n_points_x> <n_points_z> [--only-stage-1]"
  exit
fi

# directories
basedir=$(pwd)/..
opendihu_directory=$(pwd)/../../../..
stl_utility_directory=${opendihu_directory}/scripts/stl_utility
pyod=${opendihu_directory}/dependencies/python/install/bin/python3
# -----------------
# settings

# input mesh to use, this is an STL mesh of the surface of the geometry
input_file=$1

# output filenames
pickle_output_file=$2
bin_output_file=$3

# z range of the mesh for which fibers should be extracted
min_z=$4
max_z=$5

# number of fibers in x or y direction,  has to be even for quadratic elements
n_points_x=$6

# number of subdivisions of the z range, this has to be odd for quadratic elements
n_rings=$7

# choose which algorithm to use
triangulation_type=2         # triangulation_type:  0 = scipy, 1 = triangle, 2 = custom with CoG, 3 = custom with minimized distance
parametric_space_shape=3     # 0 = unit circle, 1 = unit square, 2 = unit square with adjusted grid, 3 = unit circle with adjusted grid, 4 = unit circle with optimized positions

# use the following to create it:
#$pyod ${stl_utility_directory}/remove_inside_triangles.py ${input_file} ${input_file}_repaired.stl    # remove non-surface triangles (takes ~1h)
# -----------------

# change to working directory for meshes
mkdir -p ${basedir}/processed_meshes/out
mkdir -p ${basedir}/processed_meshes/results
cd ${basedir}/processed_meshes

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

# create a mesh, reads input from `rings`, write output to `mesh` and `mesh.bin`, debugging output to out/mesh_02* - out/mesh_09* (takes ~2min)
n_grid_points_x=$(python -c "print($n_points_x+1)")

# usage: ./create_mesh.py [<triangulation_type> [<parametric_space_shape> [<n_points_x> [<n_grid_points_x> [<improve_mesh> [<pickle output filename> <bin output filename>]]]]]]
echo $pyod ${basedir}/scripts/utility/create_mesh.py $triangulation_type $parametric_space_shape $n_points_x $n_grid_points_x 1 $pickle_output_file $bin_output_file
$pyod ${basedir}/scripts/utility/create_mesh.py $triangulation_type $parametric_space_shape $n_points_x $n_grid_points_x 1 $pickle_output_file $bin_output_file

if [ "$8" = "--only-stage-1" ]; then
echo "Preprocessing complete, files ${pickle_output_file} and ${bin_output_file} were created, now exit."
exit
fi

# rename intermediate files
cp ${pickle_output_file} ../processed_meshes/results/mesh.pickle
cp ${bin_output_file} ../processed_meshes/results/mesh.bin

# ---- c++ part, streamline tracing ------
#ln -s mesh mesh_normal

# compile program
cd $opendihu_directory
scons no_tests=TRUE

# change to basedir and compile example
cd ${basedir}
scons
cd build_release

# run the simulation, created meshes with 100. elements per cm and the longest with 15 cm length
scheme_name=laplace_structured_linear

# prefix for intermediate files
prefix="created_"

# command arguments: <input_filename> <output_filename> [<target_element_length> [<target_fiber_length>]]
echo ./${scheme_name} ../settings_streamline_tracer.py ${pickle_output_file} ${prefix}${scheme_name} 1e-2 15
./${scheme_name} ../settings_streamline_tracer.py ${pickle_output_file} ${prefix}${scheme_name} 1e-2 15
exit 1

# convert csv output file to stl and pickle file
$pyod ${basedir}/scripts/utility/convert_mesh_csv_to_stl.py ${prefix}${scheme_name}.csv ../processed_meshes/results/${prefix}${scheme_name}.stl
$pyod ${basedir}/scripts/utility/convert_mesh_csv_to_pickle.py ${prefix}${scheme_name}.csv ${pickle_output_file}
$pyod ${basedir}/scripts/utility/convert_mesh_csv_to_pickle.py ${prefix}${scheme_name}_raw.csv ../processed_meshes/results/${prefix}${scheme_name}_raw.pickle

# convert pickle output file to bin fibers file
$pyod ${basedir}/scripts/utility/convert_pickle_to_bin.py ${pickle_output_file} ${bin_output_file}

