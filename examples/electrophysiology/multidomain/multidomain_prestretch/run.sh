# Create scaled geometry

export INPUT_DIR=$OPENDIHU_HOME/examples/electrophysiology/input
input_file=$INPUT_DIR/left_biceps_brachii_9x9fibers.bin
dir=$(pwd)

# get the bounding box
bounding_box_string=$(get_bounding_box_fibers.py $input_file) 
read -ra bounding_box <<< ${bounding_box_string}
zmax=${bounding_box[5]}

echo "bounding_box_string=[${bounding_box_string}], zmax=${zmax}"

# "squeeze" muscle in z direction
output_file=left_biceps_brachii_9x9fibers_shortened.bin
scale_bin_fibers.py $input_file $output_file 1 1 0.7 0 0 $zmax

# create fat mesh
fat_mesh_file=${output_file}_fat.bin
create_fat_layer.py $output_file $fat_mesh_file 0.5 3

# run simulation
example_dir=$OPENDIHU_HOME/examples/electrophysiology/multidomain/multidomain_prestretch
$example_dir/build_release/multidomain_prestretch_dummy \
  $example_dir/settings_prestretch.py $dir/shortened.py

mpirun -n 4 ./precontraction ../settings_precontraction.py precontraction.py 
mpirun -n 2 ./precontraction_prestretch ../settings_precontraction_prestretch.py precontraction_prestretch.py
mpirun -n 2 ./multidomain_prestretch ../settings_multidomain_prestretch.py multidomain.py
