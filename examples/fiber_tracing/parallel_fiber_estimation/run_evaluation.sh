for use_neumann_bc in true false; do
for improve_mesh in true false; do
for method in splines stl; do
for refinement in 1 2 3; do
  
  # clear output directory
  rm -rf out
  echo -n "$(date): " >> calls.txt 
  
  # run program
  # arguments: splines_or_stl refinement improve_mesh use_gradient_field use_neumann_bc
  echo mpirun -n 64 ./generate ../settings/settings.py $method $refinement $improve_mesh false $use_neumann_bc | tee -a calls.txt
  mpirun -n 64 ./generate ../settings/settings.py $method $refinement $improve_mesh false $use_neumann_bc
  
  # run mesh_evaluate_quality on created files
  for file in `ls -rt *.bin* | head -n 2`; do
    mesh_evaluate_quality.py $file
  done

  # rename output directory
  mv out out_${method}_${refinement}_${improve_mesh}_false_${use_neumann_bc}
  
done
done
done
done

