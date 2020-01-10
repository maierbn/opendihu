sudo perf stat -e r5301c7,r5304c7,r5310c7 > p.txt 2>&1 &
sleep 5 && sudo killall perf
compute_flops.py p.txt
#rm p.txt
