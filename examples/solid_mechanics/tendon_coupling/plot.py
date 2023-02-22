import numpy as np
import matplotlib.pyplot as plt

with open('single.txt', 'r') as f:
    lines = f.readlines()
    t = [float(line.split()[0]) for line in lines]
    z = [float(line.split()[1]) for line in lines]

with open('singlelow.txt', 'r') as f:
    lines = f.readlines()
    t_low = [float(line.split()[0]) for line in lines]
    z_low= [float(line.split()[1]) for line in lines]

with open('explicit.txt', 'r') as f:
    lines = f.readlines()
    z_explicit= [float(line.split()[1]) for line in lines]

with open('explicithigh.txt', 'r') as f:
    lines = f.readlines()
    t_high = [float(line.split()[0]) for line in lines]
    z_singlehigh= [float(line.split()[1]) for line in lines]

with open('displacement4.txt', 'r') as f:
    lines = f.readlines()
    displacement4 = [float(line.split()[1]) for line in lines]

with open('displacement6.txt', 'r') as f:
    lines = f.readlines()
    displacement6 = [float(line.split()[1]) for line in lines]

plt.figure(1)
plt.plot(t_low,z_low, "o", label="dt=0.1")
plt.plot(t,z, "o", label="dt=0.01")
plt.plot(t_high,z_singlehigh, "-", label="dt=0.001")

plt.title("Single tendon")
plt.xlabel("time (ms)")
plt.ylabel("length (cm)")
plt.legend()
plt.show()

plt.figure(2)

plt.plot(t_high,z_singlehigh, "-", label="reference")

plt.plot(t,z_explicit, "x",label="explicit")
plt.plot(t,displacement4, "o",label="implicit v1")
plt.plot(t,displacement6, "o",label="implicit v2")

plt.title("Two coupled tendons")
plt.xlabel("time (ms)")
plt.ylabel("joint length (cm)")
plt.legend()
plt.show()
