# reference
for i in `seq 1 5`; do
  aprun -n 48 -N 24 ./multiple_fibers_gnu ../settings.py 4 gnu_reference       # 4 = 4 processes per fiber
done

# packed turbo off
for i in `seq 1 5`; do
  aprun -n 48 -N 24 --p-state 2500000 ./multiple_fibers_gnu ../settings.py 4 gnu_turbo_off       # 4 = 4 processes per fiber
done

# packed slow
for i in `seq 1 5`; do
  aprun -n 48 -N 24 --p-state 1250000 ./multiple_fibers_gnu ../settings.py 4 gnu_slow       # 4 = 4 processes per fiber
done

# non-packed turbo off
for i in `seq 1 5`; do
  aprun -n 48 -N 24 -S 12 --p-state 2500000 ./multiple_fibers_gnu ../settings.py 4 gnu_non-packed_turbo_off       # 4 = 4 processes per fiber
done

# non-packed slow
for i in `seq 1 5`; do
  aprun -n 48 -N 24 -S 12 --p-state 1250000 ./multiple_fibers_gnu ../settings.py 4 gnu_non-packed_slow       # 4 = 4 processes per fiber
done

# hyperthreads
for i in `seq 1 5`; do
  aprun -n 48 -N 48 -S 24 -j 2 --p-state 1250000 ./multiple_fibers_gnu ../settings.py 4 gnu_hyperthreads       # 4 = 4 processes per fiber
done

# hyperthreads
for i in `seq 1 5`; do
  aprun -n 48 -N 48 -S 24 -j 2  ./multiple_fibers_gnu ../settings.py 4 gnu_hyperthreads_fast       # 4 = 4 processes per fiber
done


# hugepages
module unload craype-hugepages16M
module load craype-hugepages64M

for i in `seq 1 5`; do
  aprun -n 48 -N 24 ./multiple_fibers_gnu ../settings.py 4 gnu_hugepages       # 4 = 4 processes per fiber
done

module unload craype-hugepages64M
module load craype-hugepages16M
 
