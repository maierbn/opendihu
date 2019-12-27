#!/usr/bin/python3
# replace constant12] by constant12 for all numbers


input_filename = "a"
output_filename = "b"
contents = []
outfile = open(output_filename,"w")
with open(input_filename,"r") as infile:
  for line in infile:
    s = line
    search_str = "constant"
    result_line = ""
    while search_str in s:
      p = s.find(search_str)
      pend = s.find("]",p)
      ns = s[p + len(search_str):pend]
      print("s: {}".format(s))
      print("p: {}, pend: {}, ns: {}".format(p, pend, ns))
      if ns == '':
        continue
      n = (int)(ns)
      result_line += s[0:p] + "constant" + str(n)
      s = s[pend+1:]
    result_line += s
    outfile.write(result_line)
