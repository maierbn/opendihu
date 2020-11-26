#!/bin/bash

# This scripts creates the 1D fiber meshes and the 3D mesh of the muscle. 
# Input is one stl file of the geometry that was extracted from cmgui.
# Multiple files with different number of fibers are created as *.bin files that can directly be used by opendihu.
# In order to visualize a bin file, run `examine_bin_fibers.py <filename.bin>`

input_file=original_meshes/left_triceps_brachii.stl
# if you change the input file, you probably also have to experiment with the bottom_z and top_z clipping parameters for the muscle 

# [cm] range along z-axis for which the muscle volume is extracted
bottom_z_clip=2
top_z_clip=23

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
if [[ ! -f "processed_meshes/${basename}_01_no_inside_triangles.stl" ]]; then
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


# now ${basename}_03_bottom_at_zero.stl is the same as "biceps_full.stl"

echo ""
echo "--- Create a spline surface of the geometry"

# create spline surface
cd $opendihu_directory/scripts/geometry_manipulation
if [[ ! -f "${current_directory}/processed_meshes/${basename}_04_spline_surface.pickle" ]]; then
  $pyod ./create_spline_surface.py \
    ${current_directory}/processed_meshes/${basename}_03_bottom_at_zero.stl \
    ${current_directory}/processed_meshes/${basename}_04_spline_surface.stl \
    ${current_directory}/processed_meshes/${basename}_04_spline_surface.pickle \
    $bottom_z_clip $top_z_clip
else
  echo "file processed_meshes/${basename}_04_spline_surface.pickle already exists"
fi

echo ""
echo "--- Compile opendihu"
cd $opendihu_directory
$scons no_tests=TRUE

echo ""
echo "--- Compile parallel fiber estimation"
cd $parallel_fiber_estimation_directory
$scons

cd build_release

echo ""
echo "--- Generate actual fiber meshes in different sizes using the parallel_fiber_estimation example"
#read -p "Press enter to continue"

# parameters for files with different numbers of fibers
# number of fibers (number without boundary): 9(7),   11(9),   13(11),   25(23),   37(35),     65(63),     129(127),     257(255),     513(511)
# total number of fibers:                     81(49), 121(81), 169(121), 625(529), 1369(1225), 4225(3969), 16641(16129), 66049(65025), 263139(261121)
# m = number of fine grid fibers
# l = maximum recursion level
# n = number of elements in x and y coordinate directions
array_m=(0 0 0 0 0 1 1 3 7)
array_l=(0 0 0 1 1 1 2 2 2)
array_n=(4 5 6 6 9 8 8 8 8)
array_number_fibers1=(9 11 13 25 37 65 129 257 513)
array_number_fibers2=(7 9  11 23 35 63 127 255 511)

# loop over parameter combinations
for i in ${!array_l[@]}; do

  # get parameters
  l=${array_l[i]}
  m=${array_m[i]}
  n=${array_n[i]}
  number_fibers1=${array_number_fibers1[i]}
  number_fibers2=${array_number_fibers2[i]}

  # use the appropriate number of processes
  if [[ "$l" -eq "0" ]]; then
    mpi_command="mpirun -n 1"
  elif [[ "$l" -eq "1" ]]; then
    mpi_command="mpirun -n 8"
  elif [[ "$l" -eq "2" ]]; then
    mpi_command="mpirun -n 64"
  fi

  # for even number of elements use quadratic formulation
  if [[ "$n" -eq "5" || "$n" -eq "9" ]]; then
    program_name="generate"
  else
    program_name="generate_quadratic"
  fi

  echo ""
  echo "--- Generate fiber mesh files with ${number_fibers1}x${number_fibers1} and ${number_fibers2}x${number_fibers2} fibers"
  cd $parallel_fiber_estimation_directory/build_release
  
  # create file, if does not yet exist
  if [[ ! -f "${current_directory}/processed_meshes/${basename}_${number_fibers1}x${number_fibers1}fibers.bin" ]]; then
    echo "${mpi_command} ./${program_name} ../settings_generate.py \
      --input_filename_or_splines_or_stl ${current_directory}/processed_meshes/${basename}_04_spline_surface.pickle \
      --output_filename ${current_directory}/processed_meshes/${basename}_05_0x0fibers.bin \
      --bottom_z_clip $bottom_z_clip \
      --top_z_clip $top_z_clip \
      --element_size $element_length \
      --n_elements_z_per_subdomain 30 \
      --use_neumann_bc=False \
      -l=${l} -m=${m} --n_elements_x_per_subdomain=${n} \
      --program_name=${program_name}"

    ${mpi_command} ./${program_name} ../settings_generate.py \
      --input_filename_or_splines_or_stl ${current_directory}/processed_meshes/${basename}_04_spline_surface.pickle \
      --output_filename ${current_directory}/processed_meshes/${basename}_05_0x0fibers.bin \
      --bottom_z_clip $bottom_z_clip \
      --top_z_clip $top_z_clip \
      --element_size $element_length \
      --n_elements_z_per_subdomain 30 \
      --use_neumann_bc=False \
      -l=${l} -m=${m} --n_elements_x_per_subdomain=${n} \
      --program_name=${program_name}

    # move the fibers in the fibers.bin file back to their original position
    echo ""
    echo "--- Move the fibers file back to original position"
    $pyod $opendihu_directory/scripts/file_manipulation/translate_bin_fibers.py \
      ${current_directory}/processed_meshes/${basename}_05_${number_fibers1}x${number_fibers1}fibers.bin \
      ${current_directory}/processed_meshes/${basename}_06_${number_fibers1}x${number_fibers1}fibers_original_position.bin \
      0 0 ${bottom_bounding_box_value}
      
    $pyod $opendihu_directory/scripts/file_manipulation/translate_bin_fibers.py \
      ${current_directory}/processed_meshes/${basename}_05_${number_fibers2}x${number_fibers2}fibers.no_boundary.bin \
      ${current_directory}/processed_meshes/${basename}_06_${number_fibers2}x${number_fibers2}fibers_original_position.bin \
      0 0 ${bottom_bounding_box_value}

    echo ""
    echo "--- Reverse the numbering in y direction"
    $pyod $opendihu_directory/scripts/file_manipulation/reverse_x_order_bin_fibers.py \
      ${current_directory}/processed_meshes/${basename}_06_${number_fibers1}x${number_fibers1}fibers_original_position.bin \
      ${current_directory}/processed_meshes/${basename}_07_${number_fibers1}x${number_fibers1}fibers_y_reversed.bin 

    $pyod $opendihu_directory/scripts/file_manipulation/swap_xy_bin_fibers.py \
      ${current_directory}/processed_meshes/${basename}_07_${number_fibers1}x${number_fibers1}fibers_y_reversed.bin \
      ${current_directory}/processed_meshes/${basename}_08_${number_fibers1}x${number_fibers1}fibers_xy_swapped.bin 

    $pyod $opendihu_directory/scripts/file_manipulation/reverse_y_order_bin_fibers.py \
      ${current_directory}/processed_meshes/${basename}_06_${number_fibers2}x${number_fibers2}fibers_original_position.bin \
      ${current_directory}/processed_meshes/${basename}_07_${number_fibers2}x${number_fibers2}fibers_y_reversed.bin 

    $pyod $opendihu_directory/scripts/file_manipulation/swap_xy_bin_fibers.py \
      ${current_directory}/processed_meshes/${basename}_07_${number_fibers2}x${number_fibers2}fibers_y_reversed.bin \
      ${current_directory}/processed_meshes/${basename}_08_${number_fibers2}x${number_fibers2}fibers_xy_swapped.bin 

    # rename the fibers to their final name
    mv ${current_directory}/processed_meshes/${basename}_08_${number_fibers1}x${number_fibers1}fibers_xy_swapped.bin ${current_directory}/processed_meshes/${basename}_${number_fibers1}x${number_fibers1}fibers.bin
    mv ${current_directory}/processed_meshes/${basename}_08_${number_fibers2}x${number_fibers2}fibers_xy_swapped.bin ${current_directory}/processed_meshes/${basename}_${number_fibers2}x${number_fibers2}fibers.no_boundary.bin

  else
    echo "File processed_meshes/${basename}_${number_fibers1}x${number_fibers1}fibers.bin already exists, do not create again."
  fi
  
  # create fat layer meshes for the smaller files
  if [[ "${number_fibers1}" -le "100" ]]; then
  
    echo ""
    echo "--- Create fat layer meshes for ${number_fibers1}x${number_fibers1} and ${number_fibers2}x${number_fibers2} fibers"
    cd $current_directory/processed_meshes
    
    echo ""
    echo "Create ${basename}_${number_fibers1}x${number_fibers1}fibers.bin_fat.bin if does not exist"
    if [[ ! -f "${basename}_${number_fibers1}x${number_fibers1}fibers.bin_fat.bin" ]]; then
      $pyod $opendihu_directory/scripts/create_fat_layer.py ${basename}_${number_fibers1}x${number_fibers1}fibers.bin
    else
      echo "File ${basename}_${number_fibers1}x${number_fibers1}fibers.bin_fat.bin already exists, do not create again."
    fi
    
    echo ""
    echo "Create ${basename}_${number_fibers2}x${number_fibers2}fibers.no_boundary.bin_fat.bin if does not exist"
    if [[ ! -f "${basename}_${number_fibers2}x${number_fibers2}fibers.no_boundary.bin_fat.bin" ]]; then
      $pyod $opendihu_directory/scripts/create_fat_layer.py ${basename}_${number_fibers2}x${number_fibers2}fibers.no_boundary.bin
    else
      echo "File ${basename}_${number_fibers2}x${number_fibers2}fibers.no_boundary.bin_fat.bin already exists, do not create again."
    fi
  else   
    echo ""
    echo "--- Do not create a fat layer mesh for ${number_fibers1}x${number_fibers1} and ${number_fibers2}x${number_fibers2} fibers because the file would be very large."
    break
  fi

  cd $current_directory

done

cd $current_directory

