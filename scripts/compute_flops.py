#!/usr/bin/python3
#
# Parses an output file of perf and compute flops.
# Command to run perf:
#   sudo perf stat -e r5301c7,r5304c7,r5310c7 > p.txt 2>&1


import sys

filename = "p.txt"
if len(sys.argv) > 1:
  filename = sys.argv[1]

print("filename: {}".format(filename))

nOps1 = 0  # number of scalar floating point operations
nOps2 = 0  # number of sse instructions
nOps4 = 0  # number of avx-2 instructions
nSeconds = 0.0

with open(filename,"r") as f:
  for line in f:
    line = line.strip()
    if "r53" in line:
      str_number = line.split(" ")[0]
      str_number = str_number.replace(".","")
      number = (int)(str_number)
    if "r5301c7" in line:
      nOps1 = number
    elif "r5304c7" in line:
      nOps2 = number
    elif "r5310c7" in line:
      nOps2 = number

    if "seconds" in line:
      str_number = line.split(" ")[0]
      str_number = str_number.replace(",",".")
      nSeconds = (float)(str_number)

nFLOPS = (nOps1 + 2*nOps2 + 4*nOps4)/nSeconds*1e-9
print("{} Gflops, {} seconds (Skylake I5-6300U peak: 76.8 Gflops, {:.4}%)".format(nFLOPS, nSeconds, nFLOPS/76.8*100.))

