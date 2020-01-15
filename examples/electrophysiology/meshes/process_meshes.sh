#!/bin/bash

input_file=original_meshes/left_biceps_brachii.stl
# if you change the input file, you probably also have to experiment with the bottom_z and top_z clipping parameters for the muscle 

# get filename and basename
filename=${input_file##*/}    # left_biceps_brachii.stl
basename=${filename%.stl}     # left_biceps_brachii

# define directories
current_directory=$(pwd)
parallel_fiber_estimation_directory=$(pwd)/../../parallel_fiber_estimation
opendihu_directory=$(pwd)/../../..
pyod=$opendihu_directory/dependencies/python/install/bin/python3

mkdir -p processed_meshes

# remove inside triangles
echo ""
echo "--- Remove inside triangles, this will take long"
if [[ ! -f "processed_meshes/${basename}_01_no_inside_triangles.stl" ]]; then
  echo "create file \"processed_meshes/${basename}_01_no_inside_triangles.stl\""
  $pyod ./stl_utility/remove_inside_triangles.py ${input_file} processed_meshes/${basename}_01_no_inside_triangles.stl
else
  echo "File \"processed_meshes/${basename}_01_no_inside_triangles.stl\" already exists, do not create again."
fi
$pyod ./stl_utility/stl_to_binary.py \
  processed_meshes/${basename}_01_no_inside_triangles.stl \
  processed_meshes/${basename}_02_no_inside_triangles_binary.stl

# bounding box of left_biceps_brachii.stl is [36.4137001038,140.981994629] x [120.475997925,204.354003906] x [-635.0,-290.936004639]
# the z values are the negative image number from visual human male, i.e. from image 290 (shoulder) down to 635 (ellbow)
# slices have thickness of 1mm


# move mesh such that bottom is at 0
echo ""
echo "--- Move geometry such that bottom is at z=0."
bottom_bounding_box_value=`./stl_utility/get_bottom_bounding_box_value.py processed_meshes/${basename}_02_no_inside_triangles_binary.stl`
echo "Bottom bounding box value is ${bottom_bounding_box_value}"

bottom_bounding_box_value=`python -c "print(-${bottom_bounding_box_value})"`
echo "Negative bottom bounding box value is ${bottom_bounding_box_value}"  # 635.0

$pyod ./stl_utility/translate_stl.py \
  processed_meshes/${basename}_02_no_inside_triangles_binary.stl \
  processed_meshes/${basename}_03_bottom_at_zero.stl \
  0 0 ${bottom_bounding_box_value}


# now ${basename}_03_bottom_at_zero.stl is the same as "biceps_full.stl"

echo ""
echo "--- Create a spline surface of the geometry"

# create spline surface
cd $opendihu_directory/scripts/geometry_manipulation
$pyod ./create_spline_surface.py \
  ${current_directory}/processed_meshes/${basename}_03_bottom_at_zero.stl \
  ${current_directory}/processed_meshes/${basename}_04_spline_surface.stl \
  ${current_directory}/processed_meshes/${basename}_04_spline_surface.pickle \
  70 250

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
echo "--- Generate 7x7 fibers.bin file"
#read -p "Press enter to continue"

if [[ ! -f "${current_directory}/processed_meshes/${basename}_05_7x7fibers.bin" ]]; then

./generate ../settings_generate_7x7.py \
  ${current_directory}/processed_meshes/${basename}_04_spline_surface.pickle \
  ${current_directory}/processed_meshes/${basename}_05_0x0fibers.bin \
  72 220 0.1
else

echo "file processed_meshes/${basename}_05_7x7fibers.bin already exists"

fi

mv ${current_directory}/processed_meshes/7x7fibers.no_boundary.bin ${current_directory}/processed_meshes/${basename}_05_7x7fibers.bin

echo ""
echo "--- Move the fibers file back to originial position"
# move the fibers in the fibers.bin file back to their original position
$pyod $opendihu_directory/scripts/geometry_manipulation/translate_bin_fibers.py \
  ${current_directory}/processed_meshes/${basename}_05_7x7fibers.bin \
  ${current_directory}/processed_meshes/${basename}_06_7x7fibers_original_position.bin
  0 0 ${bottom_bounding_box_value}
  
cp ${current_directory}/processed_meshes/${basename}_06_7x7fibers_original_position.bin ${current_directory}/processed_meshes/${basename}_7x7fibers.bin

cd $current_directory
# refine the given, serially created file with 7x7 fibers

echo ""
echo "--- Refine fibers file to create more dense fibers"

# input fiber
input=${current_directory}/processed_meshes/${current_directory}/processed_meshes/${basename}_7x7fibers.bin

${parallel_fiber_estimation_directory}/build_release/refine ${parallel_fiber_estimation_directory}/settings_refine.py 1 $input     # 13
${parallel_fiber_estimation_directory}/build_release/refine ${parallel_fiber_estimation_directory}/settings_refine.py 3 $input     # 25
${parallel_fiber_estimation_directory}/build_release/refine ${parallel_fiber_estimation_directory}/settings_refine.py 5 $input     # 37
${parallel_fiber_estimation_directory}/build_release/refine ${parallel_fiber_estimation_directory}/settings_refine.py 10 $input     # 67
${parallel_fiber_estimation_directory}/build_release/refine ${parallel_fiber_estimation_directory}/settings_refine.py 17 $input     # 109
${parallel_fiber_estimation_directory}/build_release/refine ${parallel_fiber_estimation_directory}/settings_refine.py 30 $input     # 187
${parallel_fiber_estimation_directory}/build_release/refine ${parallel_fiber_estimation_directory}/settings_refine.py 45 $input     # 277
${parallel_fiber_estimation_directory}/build_release/refine ${parallel_fiber_estimation_directory}/settings_refine.py 70 $input     # 427
${parallel_fiber_estimation_directory}/build_release/refine ${parallel_fiber_estimation_directory}/settings_refine.py 86 $input     # 523
