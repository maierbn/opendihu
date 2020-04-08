#!/bin/bash

# This scripts creates the 1D fiber meshes and the 3D mesh of the muscle. 
# Input is one stl file of the geometry that was extracted from cmgui.
# Multiple files with different number of fibers are created as *.bin files that can directly be used by opendihu.
# In order to visualize a bin file, run `examine_bin_fibers.py <filename.bin>`

input_file=original_meshes/left_biceps_brachii.stl
# if you change the input file, you probably also have to experiment with the bottom_z and top_z clipping parameters for the muscle 

# [cm] range along z-axis for which the muscle volume is extracted
bottom_z_clip=7.2
top_z_clip=22

# [cm] length of one 1D element in z-direction, the number of elements per fiber is thus (top_z_clip-bottom_z_clip)/element_length
element_length=0.01

# get filename and basename
filename=${input_file##*/}    # left_triceps_brachii.stl
basename=${filename%.stl}     # left_triceps_brachii

# define directories
current_directory=$(pwd)
opendihu_directory=$(pwd)/../../..
parallel_fiber_estimation_directory=${opendihu_directory}/examples/fiber_tracing/parallel_fiber_estimation
stl_utility_directory=${opendihu_directory}/scripts/stl_utility
pyod=${opendihu_directory}/dependencies/python/install/bin/python3

mkdir -p processed_meshes

# bounding box of left_biceps_brachii.stl is [36.4137001038,140.981994629] x [120.475997925,204.354003906] x [-635.0,-290.936004639]
# the z values are the negative image number from visual human male, i.e. from image 290 (shoulder) down to 635 (ellbow)
# slices have thickness of 1mm

if false; then
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
if [[ ! -f "processed_meshes/${basename}_01_no_inside_triangles.stl" ]]; then
  echo "create file \"processed_meshes/${basename}_01_no_inside_triangles.stl\""
  $pyod ${stl_utility_directory}/remove_inside_triangles.py original_meshes/cm_${basename}.stl processed_meshes/${basename}_01_no_inside_triangles.stl
else
  echo "File \"processed_meshes/${basename}_01_no_inside_triangles.stl\" already exists, do not create again."
fi
$pyod ${stl_utility_directory}/stl_to_binary.py \
  processed_meshes/${basename}_01_no_inside_triangles.stl \
  processed_meshes/${basename}_02_no_inside_triangles_binary.stl

fi
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

# tendon 1
xmin=-inf
xmax=inf
ymin=-inf
ymax=inf
zmin=1.0
zmax=7.2

# extract part for first tendon
echo ""
echo "-- Extract tendon 1"
$pyod ${stl_utility_directory}/stl_extract_cuboid.py \
  processed_meshes/${basename}_03_bottom_at_zero.stl \
  processed_meshes/${basename}_04_tendon1_box.stl \
  ${xmin} ${xmax} ${ymin} ${ymax} ${zmin} ${zmax} 1


echo ""
echo "--- Create a spline surface of the geometry"

if false; then
# create spline surface
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
  ${current_directory}/processed_meshes/${basename}_05_tendon1.pickle \
  ${current_directory}/processed_meshes/${basename}_05_tendon1.bin \
  $zmin $zmax 5 10  --only-stage-1
  

echo ""
echo "--- Set seed points from muscle fibers"

muscle_fibers_file=${opendihu_directory}/examples/electrophysiology/input/left_biceps_brachii_9x9fibers.bin
  
# extract bottom nodes of fiber mesh
# arguments: set_seed_points.py <muscle_fiber_filename> <tracing_mesh_filename> <output filename> <is_bottom>
$pyod ${current_directory}/utility/set_seed_points.py \
  ${muscle_fibers_file} \
  ${current_directory}/processed_meshes/${basename}_05_tendon1.pickle \
  ${current_directory}/processed_meshes/${basename}_06_tendon1_seed_points.pickle \
  1
  
# ---- c++ part, streamline tracing ------

echo ""
echo "--- Trace streamlines"

# compile program
cd ${opendihu_directory}
scons no_tests=TRUE

# change to basedir and compile example
cd ${opendihu_directory}/examples/fiber_tracing/streamline_tracer
scons
cd build_release

# prefix for intermediate files
prefix="created_"

# command arguments: <input_filename> <output_filename> [<target_element_length> [<target_fiber_length> [<bottom_to_top>]]]
./laplace_structured_linear ../settings_streamline_tracer.py \
  ${current_directory}/processed_meshes/${basename}_06_tendon1_seed_points.pickle \
  ${current_directory}/processed_meshes/${basename}_07_tendon1_traced_fibers \
  0 0 0

# convert csv output file to stl and pickle file
$pyod ${opendihu_directory}/examples/fiber_tracing/streamline_tracer/scripts/utility/convert_mesh_csv_to_stl.py \
  ${current_directory}/processed_meshes/${basename}_07_tendon1_traced_fibers.csv \
  ${current_directory}/processed_meshes/${basename}_07_tendon1_traced_fibers.stl
  
$pyod ${opendihu_directory}/examples/fiber_tracing/streamline_tracer/scripts/utility/convert_mesh_csv_to_pickle.py \
  ${current_directory}/processed_meshes/${basename}_07_tendon1_traced_fibers.csv \
  ${current_directory}/processed_meshes/${basename}_07_tendon1_traced_fibers.pickle
  
$pyod ${opendihu_directory}/examples/fiber_tracing/streamline_tracer/scripts/utility/convert_mesh_csv_to_pickle.py \
  ${current_directory}/processed_meshes/${basename}_07_tendon1_traced_fibers_raw.csv \
  ${current_directory}/processed_meshes/${basename}_07_tendon1_traced_fibers_raw.pickle

# convert pickle output file to bin fibers file
$pyod ${opendihu_directory}/examples/fiber_tracing/streamline_tracer//scripts/utility/convert_pickle_to_bin.py 
  ${current_directory}/processed_meshes/${basename}_07_tendon1_traced_fibers.pickle
  ${current_directory}/processed_meshes/${basename}_07_tendon1_traced_fibers.bin
  
echo ""
echo "--- Move the fibers file back to original position"

# move the fibers in the fibers.bin file back to their original position
$pyod $opendihu_directory/scripts/file_manipulation/translate_bin_fibers.py \
  ${current_directory}/processed_meshes/${basename}_07_tendon1_traced_fibers.bin \
  ${current_directory}/processed_meshes/${basename}_08_tendon1_traced_fibers_original_position.bin \
  0 0 ${bottom_bounding_box_value}
  
  
cd $current_directory
  
exit
  
# ----------------------------
# tendon 2, a and b parts
xmin=-inf
xmax=inf
ymin=-inf
ymax=inf
zmin=22.9
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
echo "--- Move the file back to original position"

# move the fibers in the fibers.bin file back to their original position
$pyod ${stl_utility_directory}/translate_stl.py \
  processed_meshes/${basename}_06_tendon2a_box4.stl \
  processed_meshes/${basename}_07_tendon2a_original_position.stl \
  0 0 ${bottom_bounding_box_value}
  
echo ""
echo "--- Move the file back to original position"

# move the fibers in the fibers.bin file back to their original position
$pyod ${stl_utility_directory}/translate_stl.py \
  processed_meshes/${basename}_06_tendon2b_box4.stl \
  processed_meshes/${basename}_07_tendon2b_original_position.stl \
  0 0 ${bottom_bounding_box_value}

# do not create spline surfaces because this does not work so well for the upper tendon

echo ""
echo "--- Create a spline surface of the geometry"

if false; then

# create spline surface
cd $opendihu_directory/scripts/geometry_manipulation
$pyod ./create_spline_surface.py \
  ${current_directory}/processed_meshes/${basename}_07_tendon2a_original_position.stl \
  ${current_directory}/processed_meshes/${basename}_04_tendon2a_spline_surface.stl \
  ${current_directory}/processed_meshes/${basename}_04_tendon2a_spline_surface.pickle \
  -inf inf

# create spline surface
cd $opendihu_directory/scripts/geometry_manipulation
$pyod ./create_spline_surface.py \
  ${current_directory}/processed_meshes/${basename}_07_tendon2b_original_position.stl \
  ${current_directory}/processed_meshes/${basename}_04_tendon2b_spline_surface.stl \
  ${current_directory}/processed_meshes/${basename}_04_tendon2b_spline_surface.pickle \
  -inf inf

cd $current_directory

fi

echo ""
echo "--- Create pickle mesh"

cd $opendihu_directory/examples/fiber_tracing/streamline_tracer/scripts

# arguments <input_stl_file> <output_pickle_file> <output_bin_file> <min_z> <max_z> <n_points_x> <n_points_z> [--only_preprocessing]
./create_mesh.sh \
  ${current_directory}/processed_meshes/${basename}_04_tendon1_box.stl \
  ${current_directory}/processed_meshes/${basename}_05_tendon1.pickle \
  ${current_directory}/processed_meshes/${basename}_05_tendon1.bin \
  $zmin $zmax 5 10 
  

todo ab hier

mithilfe von streamline_tracer, run_biceps.sh mesh erzeugen, dann tracen

exit

echo ""
echo "--- Compile opendihu"
cd $opendihu_directory
scons no_tests=TRUE

echo ""
echo "--- Compile parallel fiber estimation"
cd $parallel_fiber_estimation_directory
scons

cd build_release

echo ""
echo "--- Generate 7x7 and 9x9 fibers.bin file"
#read -p "Press enter to continue"


if [[ ! -f "${current_directory}/processed_meshes/${basename}_05_7x7fibers.bin" ]]; then
  echo "./generate ../settings_generate_7x7.py ${current_directory}/processed_meshes/${basename}_04_spline_surface.pickle ${current_directory}/processed_meshes/${basename}_05_0x0fibers.bin 7.2 22 0.01"
  ./generate ../settings_generate_7x7.py \
    ${current_directory}/processed_meshes/${basename}_04_spline_surface.pickle \
    ${current_directory}/processed_meshes/${basename}_05_0x0fibers.bin \
    $bottom_z_clip $top_z_clip $element_length
else
  echo "file processed_meshes/${basename}_05_7x7fibers.bin already exists"
fi

cp ${current_directory}/processed_meshes/${basename}_05_7x7fibers.no_boundary.bin ${current_directory}/processed_meshes/${basename}_05_7x7fibers.bin

echo ""
echo "--- Move the fibers file back to original position"
# move the fibers in the fibers.bin file back to their original position
$pyod $opendihu_directory/scripts/file_manipulation/translate_bin_fibers.py \
  ${current_directory}/processed_meshes/${basename}_05_7x7fibers.bin \
  ${current_directory}/processed_meshes/${basename}_06_7x7fibers_original_position.bin \
  0 0 ${bottom_bounding_box_value}
  
$pyod $opendihu_directory/scripts/file_manipulation/translate_bin_fibers.py \
  ${current_directory}/processed_meshes/${basename}_05_9x9fibers.bin \
  ${current_directory}/processed_meshes/${basename}_06_9x9fibers_original_position.bin \
  0 0 ${bottom_bounding_box_value}
  
cp ${current_directory}/processed_meshes/${basename}_06_7x7fibers_original_position.bin ${current_directory}/processed_meshes/${basename}_7x7fibers.bin
cp ${current_directory}/processed_meshes/${basename}_06_9x9fibers_original_position.bin ${current_directory}/processed_meshes/${basename}_9x9fibers.bin

cd $current_directory
# refine the given, serially created file with 7x7 fibers

echo ""
echo "--- Refine fibers file to create more dense fibers"

# input fiber
input=${current_directory}/processed_meshes/${basename}_7x7fibers.bin

${parallel_fiber_estimation_directory}/build_release/refine ${parallel_fiber_estimation_directory}/settings_refine.py 1 $input $bottom_z_clip $top_z_clip $element_length    # 13
${parallel_fiber_estimation_directory}/build_release/refine ${parallel_fiber_estimation_directory}/settings_refine.py 3 $input $bottom_z_clip $top_z_clip $element_length     # 25
${parallel_fiber_estimation_directory}/build_release/refine ${parallel_fiber_estimation_directory}/settings_refine.py 5 $input $bottom_z_clip $top_z_clip $element_length     # 37
${parallel_fiber_estimation_directory}/build_release/refine ${parallel_fiber_estimation_directory}/settings_refine.py 10 $input $bottom_z_clip $top_z_clip $element_length     # 67
${parallel_fiber_estimation_directory}/build_release/refine ${parallel_fiber_estimation_directory}/settings_refine.py 17 $input $bottom_z_clip $top_z_clip $element_length     # 109
${parallel_fiber_estimation_directory}/build_release/refine ${parallel_fiber_estimation_directory}/settings_refine.py 30 $input $bottom_z_clip $top_z_clip $element_length     # 187
${parallel_fiber_estimation_directory}/build_release/refine ${parallel_fiber_estimation_directory}/settings_refine.py 45 $input $bottom_z_clip $top_z_clip $element_length     # 277
${parallel_fiber_estimation_directory}/build_release/refine ${parallel_fiber_estimation_directory}/settings_refine.py 70 $input $bottom_z_clip $top_z_clip $element_length     # 427
${parallel_fiber_estimation_directory}/build_release/refine ${parallel_fiber_estimation_directory}/settings_refine.py 86 $input $bottom_z_clip $top_z_clip $element_length     # 523

# create fat layer meshes
cd $current_directory/processed_meshes
$pyod $opendihu_directory/scripts/create_fat_layer.py ${basename}_7x7fibers.bin
$pyod $opendihu_directory/scripts/create_fat_layer.py ${basename}_9x9fibers.bin
$pyod $opendihu_directory/scripts/create_fat_layer.py ${basename}_13x13fibers.bin
$pyod $opendihu_directory/scripts/create_fat_layer.py ${basename}_25x25fibers.bin
$pyod $opendihu_directory/scripts/create_fat_layer.py ${basename}_37x37fibers.bin
$pyod $opendihu_directory/scripts/create_fat_layer.py ${basename}_67x67fibers.bin
$pyod $opendihu_directory/scripts/create_fat_layer.py ${basename}_109x109fibers.bin

cd $current_directory
