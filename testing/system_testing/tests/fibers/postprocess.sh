#!/bin/bash

variant="debug"
variant="release"
CHECK_RESULTS=true
mkdir -p results

# for all testcases

for name in \
  "laplace3d_structured_linear" \
  "laplace3d_structured_quadratic" \
  "laplace3d_unstructured_linear" \
  "laplace3d_unstructured_quadratic"
do

echo ""
echo "process ${name}"
echo "----------------"

# create plot "fig.pdf"
./visualize_stl.py build_${variant}/${name}.stl

# rename resulting plot and animation
mv fig.pdf results/${name}.pdf || cp ../../document/no_plot.pdf results/${name}.pdf

echo "$(date '+%d.%m.%Y %H:%M:%S'): no tests for ${name}" >> results/log_${name}.txt

# create smaller log file that only contains the 5 most recent entries
tail results/log_${name}.txt -n 5 > results/log_recent_${name}.txt

done
