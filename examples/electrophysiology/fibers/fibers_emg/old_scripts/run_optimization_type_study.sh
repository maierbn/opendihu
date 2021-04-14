# execute this script from within the build_release directory

# load modules (on pcsgs04)
module purge
module load cuda-10.2  
module load argon-tesla/gcc/10.2-openmp-custom   
module load openmpi/3.1.6-gcc-10.2 

# compile
mkorn && srr


while :; do

# run various studies
# vc
echo ===== vc 2 =====
rm -rf lib src
mpirun -n 18 ./fibers_emg ../settings_fibers_emg.py optimization_type_study.py --end_time=10 --scenario_name=vc-aovs --optimization_type=vc --use_aovs_memory_layout=true --fiber_file=../../../input/left_biceps_brachii_25x25fibers.bin
rm -rf lib src
mpirun -n 18 ./fibers_emg ../settings_fibers_emg.py optimization_type_study.py --end_time=10 --scenario_name=vc-sova --optimization_type=vc --use_aovs_memory_layout=false --fiber_file=../../../input/left_biceps_brachii_25x25fibers.bin

# vc with approximate exponential function
rm -rf lib src
mpirun -n 18 ./fibers_emg ../settings_fibers_emg.py optimization_type_study.py --end_time=10 --scenario_name=vc-aovs-apx-e --optimization_type=vc --use_aovs_memory_layout=true --approximate_exponential_function=true --fiber_file=../../../input/left_biceps_brachii_25x25fibers.bin

# simd
echo ===== simd =====
rm -rf lib src
mpirun -n 18 ./fibers_emg ../settings_fibers_emg.py optimization_type_study.py --end_time=10 --scenario_name=simd --optimization_type=simd --fiber_file=../../../input/left_biceps_brachii_25x25fibers.bin

# openmp
echo ===== openmp 6 =====
rm -rf lib src
mpirun -n 18 ./fibers_emg ../settings_fibers_emg.py optimization_type_study.py --end_time=10 --scenario_name=openmp-18-2 --optimization_type=openmp --maximum_number_of_threads=2 --fiber_file=../../../input/left_biceps_brachii_25x25fibers.bin
rm -rf lib src
mpirun -n 18 ./fibers_emg ../settings_fibers_emg.py optimization_type_study.py --end_time=10 --scenario_name=openmp-18-1 --optimization_type=openmp --maximum_number_of_threads=1 --fiber_file=../../../input/left_biceps_brachii_25x25fibers.bin
rm -rf lib src
mpirun -n 9 ./fibers_emg ../settings_fibers_emg.py optimization_type_study.py --end_time=10 --scenario_name=openmp-9-2 --optimization_type=openmp --maximum_number_of_threads=2 --fiber_file=../../../input/left_biceps_brachii_25x25fibers.bin
rm -rf lib src
mpirun -n 9 ./fibers_emg ../settings_fibers_emg.py optimization_type_study.py --end_time=10 --scenario_name=openmp-9-4 --optimization_type=openmp --maximum_number_of_threads=4 --fiber_file=../../../input/left_biceps_brachii_25x25fibers.bin
rm -rf lib src
mpirun -n 6 ./fibers_emg ../settings_fibers_emg.py optimization_type_study.py --end_time=10 --scenario_name=openmp-6-3 --optimization_type=openmp --maximum_number_of_threads=3 --fiber_file=../../../input/left_biceps_brachii_25x25fibers.bin
rm -rf lib src
mpirun -n 6 ./fibers_emg ../settings_fibers_emg.py optimization_type_study.py --end_time=10 --scenario_name=openmp-6-6 --optimization_type=openmp --maximum_number_of_threads=6 --fiber_file=../../../input/left_biceps_brachii_25x25fibers.bin

# fast-vc
echo ===== fast-vc =====
rm -rf lib src
mpirun -n 18 ./fast_fibers_emg ../settings_fibers_emg.py optimization_type_study.py --end_time=10 --scenario_name=fast-vc --optimization_type=vc --fiber_file=../../../input/left_biceps_brachii_25x25fibers.bin --enable_weak_scaling
rm -rf lib src
mpirun -n 18 ./fast_fibers_emg ../settings_fibers_emg.py optimization_type_study.py --end_time=10 --scenario_name=fast-vc-apx-e --optimization_type=vc --fiber_file=../../../input/left_biceps_brachii_25x25fibers.bin --enable_weak_scaling --approximate_exponential_function=true

# fast-gpu
echo ===== fast-gpu =====

# load gcc 11 for gpu
module purge
module load cuda-10.2 openmpi/3.1.6-gcc-10.2 argon-tesla/gcc/11-20210110-openmp
module list

rm -rf lib src
./fast_fibers_emg ../settings_fibers_emg.py optimization_type_study.py --end_time=10 --scenario_name=fast-gpu --optimization_type=gpu --fiber_file=../../../input/left_biceps_brachii_25x25fibers.bin --enable_weak_scaling

module purge
module load cuda-10.2  
module load argon-tesla/gcc/10.2-openmp-custom   
module load openmpi/3.1.6-gcc-10.2 
module list

done
