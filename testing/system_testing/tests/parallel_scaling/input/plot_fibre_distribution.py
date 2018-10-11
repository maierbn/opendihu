#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

import numpy as np
import matplotlib.pyplot as plt

file = "MU_fibre_distribution_3780.txt"
data = np.genfromtxt(file)

plt.hist(data,bins=100)
plt.xlabel("MU no")
plt.ylabel("count")

plt.savefig("histogram_fibre_distribution.pdf")
plt.show()
