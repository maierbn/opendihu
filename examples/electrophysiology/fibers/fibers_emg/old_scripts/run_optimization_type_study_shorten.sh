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
mv src/new_slow_TK_2014_12_08.c .; rm -rf lib src && mkdir src && cp new_slow_TK_2014_12_08.c src
mpirun -n 18 ./fibers_shorten_emg ../settings_fibers_emg.py shorten.py --end_time=3 --scenario_name=vc-aovs --optimization_type=vc --use_aovs_memory_layout=true   
mv src/new_slow_TK_2014_12_08.c .; rm -rf lib src && mkdir src && cp new_slow_TK_2014_12_08.c src
mpirun -n 18 ./fibers_shorten_emg ../settings_fibers_emg.py shorten.py --end_time=3 --scenario_name=vc-sova --optimization_type=vc --use_aovs_memory_layout=false   

# vc with approximate exponential function
mv src/new_slow_TK_2014_12_08.c .; rm -rf lib src && mkdir src && cp new_slow_TK_2014_12_08.c src
mpirun -n 18 ./fibers_shorten_emg ../settings_fibers_emg.py shorten.py --end_time=3 --scenario_name=vc-aovs-apx-e --optimization_type=vc --use_aovs_memory_layout=true --approximate_exponential_function=true   

# simd
echo ===== simd =====
mv src/new_slow_TK_2014_12_08.c .; rm -rf lib src && mkdir src && cp new_slow_TK_2014_12_08.c src
mpirun -n 18 ./fibers_shorten_emg ../settings_fibers_emg.py shorten.py --end_time=3 --scenario_name=simd --optimization_type=simd   

# openmp
echo ===== openmp 6 =====
mv src/new_slow_TK_2014_12_08.c .; rm -rf lib src && mkdir src && cp new_slow_TK_2014_12_08.c src
mpirun -n 18 ./fibers_shorten_emg ../settings_fibers_emg.py shorten.py --end_time=3 --scenario_name=openmp-18-2 --optimization_type=openmp --maximum_number_of_threads=2
mv src/new_slow_TK_2014_12_08.c .; rm -rf lib src && mkdir src && cp new_slow_TK_2014_12_08.c src
mpirun -n 18 ./fibers_shorten_emg ../settings_fibers_emg.py shorten.py --end_time=3 --scenario_name=openmp-18-1 --optimization_type=openmp --maximum_number_of_threads=1  
mv src/new_slow_TK_2014_12_08.c .; rm -rf lib src && mkdir src && cp new_slow_TK_2014_12_08.c src
mpirun -n 9 ./fibers_shorten_emg ../settings_fibers_emg.py shorten.py --end_time=3 --scenario_name=openmp-9-2 --optimization_type=openmp --maximum_number_of_threads=2  
mv src/new_slow_TK_2014_12_08.c .; rm -rf lib src && mkdir src && cp new_slow_TK_2014_12_08.c src
mpirun -n 9 ./fibers_shorten_emg ../settings_fibers_emg.py shorten.py --end_time=3 --scenario_name=openmp-9-4 --optimization_type=openmp --maximum_number_of_threads=4  
mv src/new_slow_TK_2014_12_08.c .; rm -rf lib src && mkdir src && cp new_slow_TK_2014_12_08.c src
mpirun -n 6 ./fibers_shorten_emg ../settings_fibers_emg.py shorten.py --end_time=3 --scenario_name=openmp-6-3 --optimization_type=openmp --maximum_number_of_threads=3  
mv src/new_slow_TK_2014_12_08.c .; rm -rf lib src && mkdir src && cp new_slow_TK_2014_12_08.c src
mpirun -n 6 ./fibers_shorten_emg ../settings_fibers_emg.py shorten.py --end_time=3 --scenario_name=openmp-6-6 --optimization_type=openmp --maximum_number_of_threads=6  

# fast-vc
echo ===== fast-vc =====
mv src/new_slow_TK_2014_12_08.c .; rm -rf lib src && mkdir src && cp new_slow_TK_2014_12_08.c src
mpirun -n 18 ./fast_fibers_shorten_emg ../settings_fibers_emg.py shorten.py --end_time=3 --scenario_name=fast-vc --optimization_type=vc  --enable_weak_scaling 
mv src/new_slow_TK_2014_12_08.c .; rm -rf lib src && mkdir src && cp new_slow_TK_2014_12_08.c src
mpirun -n 18 ./fast_fibers_shorten_emg ../settings_fibers_emg.py shorten.py --end_time=3 --scenario_name=fast-vc-apx-e --optimization_type=vc  --enable_weak_scaling --approximate_exponential_function=true

# fast-gpu
echo ===== fast-gpu =====

# load gcc 11 for gpu
module purge
module load cuda-10.2 openmpi/3.1.6-gcc-10.2 argon-tesla/gcc/11-20210110-openmp
module list

mv src/new_slow_TK_2014_12_08.c .; rm -rf lib src && mkdir src && cp new_slow_TK_2014_12_08.c src
./fast_fibers_shorten_emg ../settings_fibers_emg.py shorten.py --end_time=3 --scenario_name=fast-gpu --optimization_type=gpu  --enable_weak_scaling 

module purge
module load cuda-10.2  
module load argon-tesla/gcc/10.2-openmp-custom   
module load openmpi/3.1.6-gcc-10.2 
module list

done
