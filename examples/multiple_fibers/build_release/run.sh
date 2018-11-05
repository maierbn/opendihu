for i in `seq 1 5`; do
  aprun -n 24 ./multiple_fibers_gnu ../settings.py 4    # 4 = 4 processes per fiber
done
