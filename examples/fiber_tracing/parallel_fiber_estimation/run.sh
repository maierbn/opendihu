# settings
use_neumann_bc=true
improve_mesh=true
method=splines
refinement=1
nproc=8

# clear output directory
rm -rf out

# run program
echo mpirun -n $nproc ./generate ../settings/settings.py $method $refinement $improve_mesh false $use_neumann_bc
mpirun -n $nproc ./generate ../settings/settings.py $method $refinement $improve_mesh false $use_neumann_bc

# run mesh_evaluate_quality on created files
for file in `ls -rt *.bin* | head -n 2`; do
  mesh_evaluate_quality.py $file
done

# rename output directory
mv out out_${method}_${refinement}_${improve_mesh}_false_${use_neumann_bc}
