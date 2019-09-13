#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)
#rc('text', usetex=True)


# set global parameters for font sizes
plt.rcParams.update({'font.size': 12})
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['lines.markersize'] = 5

def fl(lambda_f):
  lambda_f_opt = 1.2
  if  .6 <= lambda_f/lambda_f_opt <= 1.4:
    return -25./4 * (lambda_f / lambda_f_opt)**2 + 25./2 * lambda_f/lambda_f_opt - 5.25
    
  return 0.0


xlist = np.linspace(0, 2.0, 100)
ylist = [fl(x) for x in xlist]

plt.plot(xlist, ylist, '-')
plt.grid(which='major')
plt.title(r'force-length relation function $f_l(\lambda_f)$')
plt.xlabel(r'$\lambda_f$ [-]')
plt.ylabel(r'force scaling [-]')

plt.savefig("fl.pdf")
plt.show()
