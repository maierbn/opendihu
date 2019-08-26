This directory contains an example of the usage of POD and reduced models.

First Step: The run schript run.sh runs the simulation to produce python outputs in the out directory of the build_debug or build_release.
Second Step: In the corresponding debug/release folder execute the check_results.py script. The snapshots are gathered in data.csv in the out directoy.
Third Step: execute the svdexample (from the corresponding debug/release folder). This reads the data.csv from the corresponding out folder and creates there the SVDresult.csv

The content of SVDresult is the transpose of the basis (transpose of the left eigenvectors)
