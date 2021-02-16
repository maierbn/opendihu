#!/bin/bash

# This scripts creates the 1D fiber meshes and the 3D mesh of the muscle. 
# Input is one stl file of the geometry that was extracted from cmgui.
# Multiple files with different number of fibers are created as *.bin files that can directly be used by opendihu.
# In order to visualize a bin file, run `examine_bin_fibers.py <filename.bin>`

input_file=original_meshes/left_biceps_brachii.stl
# if you change the input file, you probably also have to experiment with the bottom_z and top_z clipping parameters for the muscle 

# tendon 1 (bottom)
# number of points in x and z direction of the extracted mesh
tendon_1_n_points_z=11    # 11
tendon_1_n_points_x=4

# tendon 2a (top)
# number of points in x and z direction of the extracted mesh
tendon_2a_n_points_z=11   # 21
tendon_2a_n_points_x=4    # 8

# tendon 2b (top)
# number of points in x and z direction of the extracted mesh
tendon_2b_n_points_z=11   # 21
tendon_2b_n_points_x=4    # 8

# get filename and basename
filename=${input_file##*/}    # left_triceps_brachii.stl
basename=${filename%.stl}     # left_triceps_brachii

# define directories
current_directory=$(pwd)
opendihu_directory=$(pwd)/../../..
parallel_fiber_estimation_directory=${opendihu_directory}/examples/fiber_tracing/parallel_fiber_estimation
stl_utility_directory=${opendihu_directory}/scripts/stl_utility
pyod=${opendihu_directory}/dependencies/python/install/bin/python3
scons=${opendihu_directory}/dependencies/scons/scons.py

mkdir -p processed_meshes

# bounding box of left_biceps_brachii.stl is [36.4137001038,140.981994629] x [120.475997925,204.354003906] x [-635.0,-290.936004639]
# the z values are the negative image number from visual human male, i.e. from image 290 (shoulder) down to 635 (ellbow)
# slices have thickness of 1mm

# scale mesh from mm to cm, i.e. scale coordinates by factor 0.1
echo ""
echo "--- Scale mesh from mm to cm."
if [[ ! -f "original_meshes/cm_${basename}.stl" ]]; then
  echo "create file \"original_meshes/cm_${basename}.stl\""
  $pyod ${stl_utility_directory}/scale_stl.py ${input_file} original_meshes/cm_${basename}.stl 0.1 0.1 0.1
else
  echo "File \"original_meshes/cm_${basename}.stl\" already exists, do not create again."
fi

# remove inside triangles
echo ""
echo "--- Remove inside triangles, this will take long"
if [[ ! -f "processed_meshes/${basename}_01_no_inside_triangles.stl" ]] && [[ ! -f "processed_meshes/${basename}_02_no_inside_triangles_binary.stl" ]]; then
  echo "create file \"processed_meshes/${basename}_01_no_inside_triangles.stl\""
  $pyod ${stl_utility_directory}/remove_inside_triangles.py original_meshes/cm_${basename}.stl processed_meshes/${basename}_01_no_inside_triangles.stl
else
  echo "File \"processed_meshes/${basename}_01_no_inside_triangles.stl\" already exists, do not create again."
fi
$pyod ${stl_utility_directory}/stl_to_binary.py \
  processed_meshes/${basename}_01_no_inside_triangles.stl \
  processed_meshes/${basename}_02_no_inside_triangles_binary.stl

# move mesh such that bottom is at 0
echo ""
echo "--- Move geometry such that bottom is at z=0."
bottom_bounding_box_value=`${stl_utility_directory}/get_bottom_bounding_box_value.py processed_meshes/${basename}_02_no_inside_triangles_binary.stl`
echo "Bottom bounding box value is ${bottom_bounding_box_value}"

translate_value=`python -c "print(-${bottom_bounding_box_value})"`
echo "Negative bottom bounding box value is ${translate_value}"  # 63.5

$pyod ${stl_utility_directory}/translate_stl.py \
  processed_meshes/${basename}_02_no_inside_triangles_binary.stl \
  processed_meshes/${basename}_03_bottom_at_zero.stl \
  0 0 ${translate_value}

# up to here it was the same as process_meshes_biceps.sh

# ------------------------------------------------------------
# tendon 1
# bottom tendon
# ------------------------------------------------------------

xmin=-inf
xmax=inf
ymin=-inf
ymax=inf
zmin=1.0    # 4.0
zmax=7.2

# extract part for first tendon (bottom part of muscle-tendon mesh)
echo ""
echo "-- Extract tendon 1 (bottom)"
$pyod ${stl_utility_directory}/stl_extract_cuboid.py \
  processed_meshes/${basename}_03_bottom_at_zero.stl \
  processed_meshes/${basename}_04_tendon1_box.stl \
  ${xmin} ${xmax} ${ymin} ${ymax} ${zmin} ${zmax} 1


if false; then
echo ""
echo "--- Create a spline surface of the geometry"

# create spline surface - disabled
cd $opendihu_directory/scripts/geometry_manipulation
$pyod ./create_spline_surface.py \
  ${current_directory}/processed_meshes/${basename}_04_tendon1_box.stl \
  ${current_directory}/processed_meshes/${basename}_04_tendon1_spline_surface.stl \
  ${current_directory}/processed_meshes/${basename}_04_tendon1_spline_surface.pickle \
  ${zmin} ${zmax}
fi
cd $current_directory

echo ""
echo "--- Create pickle mesh"

cd $opendihu_directory/examples/fiber_tracing/streamline_tracer/scripts

# arguments <input_stl_file> <output_pickle_file> <output_bin_file> <min_z> <max_z> <n_points_x> <n_points_z> [--only-stage-1]
./create_mesh.sh \
  ${current_directory}/processed_meshes/${basename}_04_tendon1_box.stl \
  ${current_directory}/processed_meshes/${basename}_05_tendon1_9x9.pickle \
  ${current_directory}/processed_meshes/${basename}_05_tendon1_9x9.bin \
  $zmin $zmax $tendon_1_n_points_x $tendon_1_n_points_z --only-stage-1

# transform the bin file to a vts file for debugging
echo ""
echo "--- for debugging, create vts file which can be viewed by ParaView"
cp ${current_directory}/processed_meshes/${basename}_05_tendon1_9x9.bin ${current_directory}/processed_meshes/${basename}_06_tendon1_9x9.bin
$pyod ${opendihu_directory}/scripts/file_manipulation/examine_bin_fibers.py \
  ${current_directory}/processed_meshes/${basename}_06_tendon1_9x9.bin


echo ""
echo "--- move tendon file back to original position"
  
# move tendon file back to original position
$pyod ${opendihu_directory}/scripts/file_manipulation/translate_bin_fibers.py \
  ${current_directory}/processed_meshes/${basename}_06_tendon1_9x9.bin \
  ${current_directory}/processed_meshes/${basename}_07_tendon1_9x9_original_position.bin \
  0 0 ${bottom_bounding_box_value}

echo ""
echo "--- Reverse the numbering in y direction"
$pyod ${opendihu_directory}/scripts/file_manipulation/reverse_y_order_bin_fibers.py \
  ${current_directory}/processed_meshes/${basename}_07_tendon1_9x9_original_position.bin \
  ${current_directory}/processed_meshes/${basename}_08_tendon1_9x9_y_reversed.bin

echo ""
echo "--- Swap x and y coordinates in numbering"
$pyod ${opendihu_directory}/scripts/file_manipulation/swap_xy_bin_fibers.py \
  ${current_directory}/processed_meshes/${basename}_08_tendon1_9x9_y_reversed.bin \
  ${current_directory}/processed_meshes/${basename}_08_tendon1_9x9_xy_swapped.bin

if false; then   # disable because fiber tracing is not used up to now
echo ""
echo "--- adjust top layer of nodes of tendon to match muscle file"
muscle_fibers_file=${opendihu_directory}/examples/electrophysiology/input/left_biceps_brachii_9x9fibers.bin

# set_nodes_to_match_other_mesh.py <input_filename mesh where to change nodes> <input_filename mesh to which to adapt> <output filename> <is_bottom if change bottom layer of mesh else top layer>
$pyod ${current_directory}/utility/set_nodes_to_match_other_mesh.py \
  ${current_directory}/processed_meshes/${basename}_08_tendon1_9x9_xy_swapped.bin \
  ${muscle_fibers_file} \
  ${current_directory}/processed_meshes/${basename}_09_tendon1_9x9.bin\
  0
  
fi

# transform the bin file to a vts file for debugging
echo ""
echo "--- for debugging, create vts file which can be viewed by ParaView"
$pyod ${opendihu_directory}/scripts/file_manipulation/examine_bin_fibers.py \
  ${current_directory}/processed_meshes/${basename}_08_tendon1_9x9_xy_swapped.bin
#$pyod ${opendihu_directory}/scripts/file_manipulation/examine_bin_fibers.py \
#  ${current_directory}/processed_meshes/${basename}_09_tendon1_9x9.bin
  
if false; then   # disabled because we don't do fiber tracing for now
# ---- c++ part, streamline tracing ------

echo ""
echo "--- Trace streamlines"

# compile program
cd ${opendihu_directory}
$scons no_tests=TRUE

# change to basedir and compile example
cd ${opendihu_directory}/examples/fiber_tracing/streamline_tracer
$scons
cd build_release

# prefix for intermediate files
prefix="created_"

# command arguments: <input_filename> <output_filename> [<target_element_length> [<target_fiber_length> [<bottom_to_top>]]]
./laplace_structured_linear ../settings_streamline_tracer.py \
  ${current_directory}/processed_meshes/${basename}_08_tendon1_seed_points.pickle \
  ${current_directory}/processed_meshes/${basename}_09_tendon1_traced_fibers \
  0 0 0

# convert csv output file to stl and pickle file
$pyod ${opendihu_directory}/examples/fiber_tracing/streamline_tracer/scripts/utility/convert_mesh_csv_to_stl.py \
  ${current_directory}/processed_meshes/${basename}_09_tendon1_traced_fibers.csv \
  ${current_directory}/processed_meshes/${basename}_09_tendon1_traced_fibers.stl
  
$pyod ${opendihu_directory}/examples/fiber_tracing/streamline_tracer/scripts/utility/convert_mesh_csv_to_pickle.py \
  ${current_directory}/processed_meshes/${basename}_09_tendon1_traced_fibers.csv \
  ${current_directory}/processed_meshes/${basename}_09_tendon1_traced_fibers.pickle
  
$pyod ${opendihu_directory}/examples/fiber_tracing/streamline_tracer/scripts/utility/convert_mesh_csv_to_pickle.py \
  ${current_directory}/processed_meshes/${basename}_09_tendon1_traced_fibers_raw.csv \
  ${current_directory}/processed_meshes/${basename}_09_tendon1_traced_fibers_raw.pickle

# convert pickle output file to bin fibers file
$pyod ${opendihu_directory}/examples/fiber_tracing/streamline_tracer//scripts/utility/convert_pickle_to_bin.py 
  ${current_directory}/processed_meshes/${basename}_09_tendon1_traced_fibers.pickle
  ${current_directory}/processed_meshes/${basename}_09_tendon1_traced_fibers.bin
  
echo ""
echo "--- Move the fibers file back to original position"

# move the fibers in the fibers.bin file back to their original position
$pyod $opendihu_directory/scripts/file_manipulation/translate_bin_fibers.py \
  ${current_directory}/processed_meshes/${basename}_09_tendon1_traced_fibers.bin \
  ${current_directory}/processed_meshes/${basename}_10_tendon1_traced_fibers_original_position.bin \
  0 0 ${bottom_bounding_box_value}
  
fi  # endif streamline tracing

# copy resulting mesh to input folder
cp ${current_directory}/processed_meshes/${basename}_08_tendon1_9x9_y_reversed.bin ${current_directory}/processed_meshes/${basename}_tendon1.bin
  
cd $current_directory

# ------------------------------------------------------------
# tendons 2a and 2b
# top tendons
# ------------------------------------------------------------
xmin=-inf
xmax=inf
ymin=-inf
ymax=inf
zmin=22
zmax=30

# extract part for tendon 2
echo ""
echo "-- Extract tendon 2"
$pyod ${stl_utility_directory}/stl_extract_cuboid.py \
  processed_meshes/${basename}_03_bottom_at_zero.stl \
  processed_meshes/${basename}_04_tendon2_box1.stl \
  ${xmin} ${xmax} ${ymin} ${ymax} ${zmin} ${zmax} 1

rotation_point="11.88942022 18.36920895 29.51913578"

echo ""
echo "-- rotate"
$pyod ${stl_utility_directory}/rotate_stl.py \
  processed_meshes/${basename}_04_tendon2_box1.stl \
  processed_meshes/${basename}_05_tendon2_box2.stl \
  0 -18 0 $rotation_point

xsplit=12.23
echo ""
echo "-- extract part1"
$pyod ${stl_utility_directory}/stl_extract_cuboid.py \
  processed_meshes/${basename}_05_tendon2_box2.stl \
  processed_meshes/${basename}_06_tendon2a_box3.stl \
  -inf $xsplit -inf inf -inf inf 1

echo ""
echo "-- extract part2"
$pyod ${stl_utility_directory}/stl_extract_cuboid.py \
  processed_meshes/${basename}_05_tendon2_box2.stl \
  processed_meshes/${basename}_06_tendon2b_box3.stl \
  $xsplit inf -inf inf -inf inf 2

echo ""
echo "-- rotate back"
$pyod ${stl_utility_directory}/rotate_stl.py \
  processed_meshes/${basename}_06_tendon2a_box3.stl \
  processed_meshes/${basename}_06_tendon2a_box4.stl \
  0 18 0 $rotation_point
  
echo ""
echo "-- rotate back"
$pyod ${stl_utility_directory}/rotate_stl.py \
  processed_meshes/${basename}_06_tendon2b_box3.stl \
  processed_meshes/${basename}_06_tendon2b_box4.stl \
  0 18 0 $rotation_point

echo ""
echo "--- Create pickle mesh"

cd $opendihu_directory/examples/fiber_tracing/streamline_tracer/scripts

# arguments <input_stl_file> <output_pickle_file> <output_bin_file> <min_z> <max_z> <n_points_x> <n_points_z> [--only-stage-1]
echo " --- for tendon2a"
./create_mesh.sh \
  ${current_directory}/processed_meshes/${basename}_06_tendon2a_box4.stl \
  ${current_directory}/processed_meshes/${basename}_07_tendon2a_9x9.pickle \
  ${current_directory}/processed_meshes/${basename}_07_tendon2a_9x9.bin \
  $zmin $zmax $tendon_2a_n_points_x $tendon_2a_n_points_z --only-stage-1

echo " --- for tendon2b"  
./create_mesh.sh \
  ${current_directory}/processed_meshes/${basename}_06_tendon2b_box4.stl \
  ${current_directory}/processed_meshes/${basename}_07_tendon2b_9x9.pickle \
  ${current_directory}/processed_meshes/${basename}_07_tendon2b_9x9.bin \
  $zmin $zmax $tendon_2b_n_points_x $tendon_2b_n_points_z --only-stage-1

# transform the bin file to a vts file for debugging
echo ""
echo "--- for debugging, create vts file which can be viewed by ParaView"
cp ${current_directory}/processed_meshes/${basename}_07_tendon2a_9x9.bin ${current_directory}/processed_meshes/${basename}_08_tendon2a_9x9.bin
$pyod ${opendihu_directory}/scripts/file_manipulation/examine_bin_fibers.py \
  ${current_directory}/processed_meshes/${basename}_08_tendon2a_9x9.bin
cp ${current_directory}/processed_meshes/${basename}_07_tendon2b_9x9.bin ${current_directory}/processed_meshes/${basename}_08_tendon2b_9x9.bin
$pyod ${opendihu_directory}/scripts/file_manipulation/examine_bin_fibers.py \
  ${current_directory}/processed_meshes/${basename}_08_tendon2b_9x9.bin

echo ""
echo "--- move tendon file back to original position"
  
# move tendon file back to original position
$pyod ${opendihu_directory}/scripts/file_manipulation/translate_bin_fibers.py \
  ${current_directory}/processed_meshes/${basename}_08_tendon2a_9x9.bin \
  ${current_directory}/processed_meshes/${basename}_09_tendon2a_9x9_original_position.bin \
  0 0 ${bottom_bounding_box_value}

# move tendon file back to original position
$pyod ${opendihu_directory}/scripts/file_manipulation/translate_bin_fibers.py \
  ${current_directory}/processed_meshes/${basename}_08_tendon2b_9x9.bin \
  ${current_directory}/processed_meshes/${basename}_09_tendon2b_9x9_original_position.bin \
  0 0 ${bottom_bounding_box_value}

echo ""
echo "--- Reverse the numbering in y direction"
$pyod ${opendihu_directory}/scripts/file_manipulation/reverse_y_order_bin_fibers.py \
  ${current_directory}/processed_meshes/${basename}_09_tendon2a_9x9_original_position.bin \
  ${current_directory}/processed_meshes/${basename}_10_tendon2a_y_reversed.bin
  
$pyod ${opendihu_directory}/scripts/file_manipulation/swap_xy_bin_fibers.py \
  ${current_directory}/processed_meshes/${basename}_10_tendon2a_y_reversed.bin \
  ${current_directory}/processed_meshes/${basename}_10_tendon2a_xy_swapped.bin
  
cp ${current_directory}/processed_meshes/${basename}_10_tendon2a_xy_swapped.bin ${current_directory}/processed_meshes/${basename}_tendon2a.bin

$pyod ${opendihu_directory}/scripts/file_manipulation/reverse_y_order_bin_fibers.py \
  ${current_directory}/processed_meshes/${basename}_09_tendon2b_9x9_original_position.bin \
  ${current_directory}/processed_meshes/${basename}_10_tendon2b_y_reversed.bin
  
$pyod ${opendihu_directory}/scripts/file_manipulation/swap_xy_bin_fibers.py \
  ${current_directory}/processed_meshes/${basename}_10_tendon2b_y_reversed.bin \
  ${current_directory}/processed_meshes/${basename}_10_tendon2b_xy_swapped.bin
  
cp ${current_directory}/processed_meshes/${basename}_10_tendon2b_xy_swapped.bin ${current_directory}/processed_meshes/${basename}_tendon2b.bin

cd $current_directory


# copy resulting meshes to input folder
cp ${current_directory}/processed_meshes/${basename}_tendon1.bin ${opendihu_directory}/examples/electrophysiology/input
cp ${current_directory}/processed_meshes/${basename}_tendon2a.bin ${opendihu_directory}/examples/electrophysiology/input
cp ${current_directory}/processed_meshes/${basename}_tendon2b.bin ${opendihu_directory}/examples/electrophysiology/input
