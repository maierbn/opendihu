for i in `seq 1 5`; do
  aprun -n 48 -N 24 ./multiple_fibers_gnu ../settings.py 4 gnu_reference       # 4 = 4 processes per fiber
done
