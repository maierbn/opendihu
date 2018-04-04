#!/bin/bash

variant="debug"

# for all testcases

#for name in "2d_structured_regular_fixed_linear"
#for name in \
#  "1d_structured_regular_fixed_linear" \
#  "1d_structured_regular_fixed_quadratic" \
#  "1d_structured_regular_fixed_hermite" \
#  "1d_structured_deformable_linear" \
#  "1d_structured_deformable_linear" \
#  "1d_structured_deformable_hermite" \
#  "1d_unstructured_deformable_linear" \
#  "1d_unstructured_deformable_quadratic" \
#  "1d_unstructured_deformable_hermite" \
#  \
#  "2d_structured_regular_fixed_linear" \
#  "2d_structured_regular_fixed_quadratic" \
#  "2d_structured_regular_fixed_hermite" \
#  "2d_structured_deformable_linear" \
#  "2d_structured_deformable_linear" \
#  "2d_structured_deformable_hermite" \
#  "2d_unstructured_deformable_linear" \
#  "2d_unstructured_deformable_quadratic" \
#  "2d_unstructured_deformable_hermite" \
#  \
#  "3d_structured_regular_fixed_linear" \
#  "3d_structured_regular_fixed_quadratic" \
#  "3d_structured_regular_fixed_hermite" \
#  "3d_structured_deformable_linear" \
#  "3d_structured_deformable_linear" \
#  "3d_structured_deformable_hermite" \
#  "3d_unstructured_deformable_linear" \
#  "3d_unstructured_deformable_quadratic" \
#  "3d_unstructured_deformable_hermite"
for name in \
  "2d_structured_regular_fixed_linear" \
  "2d_structured_regular_fixed_quadratic" \
  "2d_structured_regular_fixed_hermite" \
  "2d_structured_deformable_linear" \
  "2d_structured_deformable_linear" \
  "2d_structured_deformable_hermite" \
  "2d_unstructured_deformable_linear" \
  "2d_unstructured_deformable_quadratic" \
  "2d_unstructured_deformable_hermite"
do

# arguments to plot.py and check_results.py: <1=show plot window, 0=don't> <filenames>

# create plots
#../../../../scripts/plot.py 1 build_${variant}/out/${name}*

# compare to analytical solution and check if tests pass or fail (also creates an animation file "numerical_analytical.mp4")
python check_results.py 0 build_${variant}/out/${name}*

# move resulting animation
mkdir -p results
mv numerical_analytical.mp4 results/${name}_numerical_analytical.mp4
mv log.txt results/log.txt

done

#for cmguifile in `ls results/*.com`
#do
#    # run cmgui script file and export screenshot
#    echo $cmguifile
#    cmgui $cmguifile
#done
