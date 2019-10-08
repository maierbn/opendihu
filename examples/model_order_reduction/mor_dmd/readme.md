This directory contains an example of Model Order Reduction utilizing Dynamic Mode Decomposition.

First Step: The run schript run.sh runs the simulation to produce python outputs in the out directory of the build_debug or build_release.
Second Step: In the corresponding debug/release folder execute the check_results.py script. The snapshots are gathered in data.csv in the out directoy.
Third Step: execute the DMDexample (from the corresponding debug/release folder). This reads the data.csv from the corresponding out folder and creates there the DMDresult.csv

The content of DMDresult is the 
