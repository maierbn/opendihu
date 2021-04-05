export OPENDIHU_HOME=/data/scratch/sgs/maierbn/opendihu


# generate fat mesh
# usage: ./create_fat_layer.py [<fibers input file> [<output file of fat tissue> [<thickness in cm> [<y_size>]]]]
#$OPENDIHU_HOME/dependencies/python/install/bin/python3 \
#   $OPENDIHU_HOME/scripts/create_fat_layer.py \
#   $OPENDIHU_HOME/examples/electrophysiology/input/left_biceps_brachii_37x37fibers.bin \
#   $OPENDIHU_HOME/examples/electrophysiology/input/left_biceps_brachii_37x37fibers_thin_fat.bin \
#   0.2 3

# generate fiber distribution
#generate_fiber_distribution.py MU_fibre_distribution_37x37_20a.txt 20 1 37 1.05

# run
export EXAMPLE_DIR=$OPENDIHU_HOME/examples/electrophysiology/fibers/fibers_fat_emg
export dir=$(pwd)

mpirun -n 16 $EXAMPLE_DIR/build_release/fibers_fat_emg $EXAMPLE_DIR/settings_fibers_fat_emg.py $dir/20mus_better_emg.py

