#!/bin/bash

variant="debug"
CHECK_RESULTS=false
mkdir -p results

# for all testcases

#for name in "1d_structured_regular_fixed_linear"
for name in \
  "1d_structured_regular_fixed_linear" \
  "1d_structured_regular_fixed_quadratic" \
  "1d_structured_regular_fixed_hermite" \
  "1d_structured_deformable_linear" \
  "1d_structured_deformable_quadratic" \
  "1d_structured_deformable_hermite" \
  "1d_unstructured_deformable_linear" \
  "1d_unstructured_deformable_quadratic" \
  "1d_unstructured_deformable_hermite" \
  \
  "2d_structured_regular_fixed_linear" \
  "2d_structured_regular_fixed_quadratic" \
  "2d_structured_regular_fixed_hermite" \
  "2d_structured_deformable_linear" \
  "2d_structured_deformable_quadratic" \
  "2d_structured_deformable_hermite" \
  "2d_unstructured_deformable_linear" \
  "2d_unstructured_deformable_quadratic" \
  "2d_unstructured_deformable_hermite" \
  \
  "3d_structured_regular_fixed_linear" \
  "3d_structured_regular_fixed_quadratic" \
  "3d_structured_regular_fixed_hermite" \
  "3d_structured_deformable_linear" \
  "3d_structured_deformable_quadratic" \
  "3d_structured_deformable_hermite" \
  "3d_unstructured_deformable_linear" \
  "3d_unstructured_deformable_quadratic" \
  "3d_unstructured_deformable_hermite"
#for name in \
#  "1d_structured_regular_fixed_linear" \
#  "1d_structured_regular_fixed_quadratic" \
#  "1d_structured_regular_fixed_hermite" \
#  "1d_structured_deformable_linear" \
#  "1d_structured_deformable_linear" \
#  "1d_structured_deformable_hermite" \
#  "1d_unstructured_deformable_linear" \
#  "1d_unstructured_deformable_quadratic" \
#  "1d_unstructured_deformable_hermite"
do

echo "check ${name}"

# arguments to plot.py and check_results.py: <1=show plot window, 0=don't> <filenames>

# create animation "anim.mp4" and plot "fig.pdf"
../../../../scripts/plot.py 0 build_${variant}/out/${name}*

# rename resulting plot and animation
mv fig.pdf results/${name}.pdf || cp ../../document/no_plot.pdf results/${name}.pdf
mv anim.mp4 results/${name}.mp4

if [ "$CHECK_RESULTS" = true ] ; then
    
  # compare to analytical solution and check if tests pass or fail (also creates an animation file "numerical_analytical.mp4")
  mv results/log_${name}.txt log.txt
  python check_results.py 0 build_${variant}/out/${name}*
  mv log.txt results/log_${name}.txt

  # Note the behaviour of check_results.py with respect to log files:
  # It appends messages to the end of a file log.txt. 
  # Therefore we mv the named log file for the testcase to log.txt beforehand,
  # let check_results.py append to it and rename it back afterwards.

  # move resulting animation
  mv numerical_analytical.mp4 results/${name}_numerical_analytical.mp4
else
  
  echo "$(date '+%d.%m.%Y %H:%M:%S'): test ${name} disabled" >> results/log_${name}.txt
  
fi

# create smaller log file that only contains the 5 most recent entries
tail results/log_${name}.txt -n 5 > results/log_recent_${name}.txt

done
