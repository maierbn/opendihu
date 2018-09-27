#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

import numpy as np
import matplotlib.pyplot as plt

firing_times_file = "../input/MU_firing_times_real.txt"
firing_times = np.genfromtxt(firing_times_file)

print firing_times
data = np.sum(firing_times,0)
print data

plt.bar(range(len(data)),data,align='center')
plt.xlabel("MU no")
plt.ylabel("count")

plt.savefig("histogram_firing_times.pdf")
plt.show()
