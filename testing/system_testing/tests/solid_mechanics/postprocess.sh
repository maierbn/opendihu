#!/bin/bash

# set to debug or release
variant="debug"
#variant="release"
mkdir -p results


  #"mooney_rivlin_incompressible_penalty2d_numeric_jacobian_scenario_1" \
  #"mooney_rivlin_incompressible_penalty2d_analytic_jacobian_scenario_1" \
  #"mooney_rivlin_incompressible_penalty2d_numeric_jacobian_scenario_2" \
  #"mooney_rivlin_incompressible_penalty2d_analytic_jacobian_scenario_2" \
  #"mooney_rivlin_incompressible_penalty2d_numeric_jacobian_scenario_3" \
  #"mooney_rivlin_incompressible_penalty2d_analytic_jacobian_scenario_3"
  #"mooney_rivlin_incompressible_mixed_condensation" \
  #"mooney_rivlin_incompressible_mixed"

# for all testcases
for name in \
  "mooney_rivlin_incompressible_mixed2d_numeric_jacobian_scenario_2" \
  "mooney_rivlin_incompressible_mixed2d_analytic_jacobian_scenario_2" \
do

echo "check ${name}"

# arguments to plot.py and check_results.py: <1=show plot window, 0=don't> <filenames>

# create plot with name "fig.pdf"
../../../../scripts/plot.py 0 build_${variant}/out/${name}*

# rename resulting plot
mv fig.pdf results/${name}.pdf

# compare to analytical solution and check if tests pass or fail
mv results/log_${name}.txt log.txt
python check_results.py 0 build_${variant}/out/${name}*
mv log.txt results/log_${name}.txt

# create smaller log file that only contains the 5 most recent entries
tail results/log_${name}.txt -n 5 > results/log_recent_${name}.txt

# Note the behaviour of check_results.py with respect to log files:
# It appends messages to the end of a file log.txt. 
# Therefore we mv the named log file for the testcase to log.txt beforehand,
# let check_results.py append to it and rename it back afterwards.


done  # end of for loop
