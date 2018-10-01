#!/bin/bash

variant="debug"
variant="release"
CHECK_RESULTS=true
mkdir -p results

echo ""
echo "postprocess monodomain_timestep_widths"
echo "----------------"

# create animation "anim.mp4" and plot "fig.pdf"
plot.py build_${variant}

# rename resulting plot and animation
mv fig.pdf results/${name}.pdf || cp ../../document/no_plot.pdf results/${name}.pdf
mv anim.mp4 results/${name}.mp4

echo "$(date '+%d.%m.%Y %H:%M:%S'): no tests for ${name}" >> results/log_${name}.txt

# create smaller log file that only contains the 5 most recent entries
tail results/log_${name}.txt -n 5 > results/log_recent_${name}.txt
