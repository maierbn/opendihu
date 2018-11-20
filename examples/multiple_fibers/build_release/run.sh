for i in `seq 1 5`; do
  aprun -n 48 ./multiple_fibers ../settings.py 12    # 4 = 4 processes per fiber
done
