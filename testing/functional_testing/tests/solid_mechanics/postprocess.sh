#!/bin/bash

variant="debug"

# for all testcases
for name in \
  "mooney_rivlin_incompressible_mixed_condensation" \
  "mooney_rivlin_incompressible_penalty" \
  "mooney_rivlin_incompressible_mixed"
do

# arguments to plot.py and check_results.py: <1=show plot window, 0=don't> <filenames>

# create plots
#../../../../scripts/plot.py 1 build_${variant}/out/${name}*

# compare to analytical solution and check if tests pass or fail (also creates an animation file "numerical_analytical.mp4")
python check_results.py 0 build_${variant}/out/${name}*

# move resulting animation
mkdir -p results
mv numerical_analytical.mp4 results/${name}_numerical_analytical.mp4

done

#for cmguifile in `ls results/*.com`
#do
#    # run cmgui script file and export screenshot
#    echo $cmguifile
#    cmgui $cmguifile
#done
