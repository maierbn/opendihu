#!/usr/bin/python3
import numpy as np

def log(x_value):
  x = x_value-1
  x2 = x*x
  a0 = -np.log(2); a1 = 2;a2 = -1;a3 = 2/3.;a4 = -0.5;
  T0 = 1; T1 = x; T2 = 2*x2 - 1; T3 = 4*x2*x - 3*x; T4 = 8*x2*x2 - 8*x2 + 1;
  return a0*T0 + a1*T1 + a2*T2 + a3*T3 + a4*T4

def log2(x_value):
  if x_value < 2:
    t = x_value-1;
    t2 = t*t
    t4 = t2*t2
    t6 = t4*t2
    return  t - 0.5*t2 + 1./3*t2*t - 0.25*t4 + 0.2*t4*t - 1./6*t6
  elif x_value < 6:
    t = x_value-3
    t2 = t*t
    t4 = t2*t2
    t6 = t4*t2
    return  np.log(3) + 1./3*t - 1./18*t2 + 1./81*t2*t - 1/324.*t4 + 1./1215*t4*t - 1./4374*t6
  else:
    t = x_value-9
    t2 = t*t
    t4 = t2*t2
    t6 = t4*t2
    return  np.log(9) + 1./9*t - 1./162*t2 + 1./2187*t2*t - 1/26244.*t4 + 1./295245*t4*t - 1./3188646*t6

import matplotlib.pyplot as plt
fig = plt.figure()
x_list = np.linspace(0.2, 20, 200)

plt.plot(x_list, np.log(x_list), label="exact")
#plt.plot(x_list, [log(x) for x in x_list], label="Chebyshev")
plt.plot(x_list, [log2(x) for x in x_list], label="Taylor")
plt.legend()

# define global plotting parameters
plt.rcParams.update({'font.size': 16})
plt.rcParams['lines.linewidth'] = 3

fig = plt.figure()
plt.plot(x_list, [(np.log(x)-log2(x))/np.log(x) for x in x_list])
plt.grid()
plt.savefig("apxlog.pdf")
plt.show()

max_rel_error = 0
location = 0
for x in np.linspace(0.2,20,200):
  rel_error =  (np.log(x)-log2(x))/np.log(x)
  print("{}, {}, {}, error: {}".format(x, np.log(x), log(x), rel_error))
  old_rel_error = max_rel_error
  max_rel_error = max(max_rel_error, rel_error)
  if old_rel_error != max_rel_error:
    location = x
  
print("maximum relative error: {} at {}".format(max_rel_error, location))

plt.show()
