# refine the given, serially created file with 7x7 fibers

# input fiber
#input=/store/software/opendihu/examples/electrophysiology/input/7x7fibers.bin
input=/data/scratch/maierbn/input/7x7fibers.bin

cd build_release
./refine ../settings_refine.py 1 $input     # 13
./refine ../settings_refine.py 3 $input     # 25
./refine ../settings_refine.py 5 $input     # 37
./refine ../settings_refine.py 10 $input     # 67
./refine ../settings_refine.py 17 $input     # 109
./refine ../settings_refine.py 30 $input     # 187
./refine ../settings_refine.py 45 $input     # 277
./refine ../settings_refine.py 70 $input     # 427
./refine ../settings_refine.py 86 $input     # 523
cd -
