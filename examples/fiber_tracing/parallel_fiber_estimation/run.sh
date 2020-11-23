# This script runs the (potentially parallel) algorithm to generate a 3D mesh and 1D fiber meshes for the biceps muscle.
# Input is a representation of the biceps surface which is already given in the input directory.
# Output is a *.bin file that contains the points of the fibers.

# parameters for this script
use_neumann_bc=false      # (false,true) set to false to use Dirichlet BC (which is better)
improve_mesh=true         # (false,true) if the Laplacian smoothing and fixing of invalid quadrilaterals should be enabled
method=splines            # (splines,stl) which input to use, either splines = use "biceps_splines.stl" or stl = use "biceps.surface.pickle"
refinement=1              # (0,1,2,...) how often to refine the mesh for the Laplace problem, 1=twice as many elements per direction, 2=4x as many elements per direction
nproc=1                   # (1,8,64,...) number of processes, should be 8^i
m=2                       # (0,1,2,...) parameter m of fine grid fibers, a higher value means more fibers, use "compute_sizes.py" to get an estimate

# clear output directory
rm -rf out

# run program with MPI
mpirun -n $nproc --oversubscribe ./generate ../settings_generate.py \
  --input_filename_or_splines_or_stl $method \
  --refinement_factor $refinement \
  --improve_mesh $improve_mesh \
  --use_neumann_bc $use_neumann_bc \
  -m $m

# run mesh_evaluate_quality on created files
for file in `ls -rt *.bin | head -n 2`; do
  echo ""
  mesh_evaluate_quality.py $file
done

# rename output directory
mv out out_${method}_${refinement}_${improve_mesh}_${use_neumann_bc}_m${m}
