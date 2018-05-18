#!/bin/bash

variant="debug"
#variant="release"
CHECK_RESULTS=true
mkdir -p results

# create animation "anim.mp4" and plot "fig.pdf"
../../../../scripts/plot.py 0 build_${variant}/out/${name}*

# rename resulting plot and animation
mv fig.pdf results/${name}.pdf || cp ../../document/no_plot.pdf results/${name}.pdf
mv anim.mp4 results/${name}.mp4

  
echo "$(date '+%d.%m.%Y %H:%M:%S'): test ${name} disabled" >> results/log_${name}.txt
  
done
