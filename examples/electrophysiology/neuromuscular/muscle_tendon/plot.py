import numpy as np
import matplotlib.pyplot as plt

with open('single_muscle.txt', 'r') as f:
    lines = f.readlines()
    t = [float(line.split()[0]) for line in lines]
    z = [float(line.split()[1]) for line in lines]


with open('single_muscle_low.txt', 'r') as f:
    lines = f.readlines()
    t_low = [float(line.split()[0]) for line in lines]
    z_low = [float(line.split()[1]) for line in lines]

with open('single_muscle_high.txt', 'r') as f:
    lines = f.readlines()
    t_high = [float(line.split()[0]) for line in lines]
    z_high = [float(line.split()[1]) for line in lines]

with open('muscle_with_large_tendon_explicit.txt', 'r') as f:
    lines = f.readlines()
    t_large= [float(line.split()[0]) for line in lines]
    z_large= [float(line.split()[1]) for line in lines]

with open('muscle_with_small_tendon_explicit.txt', 'r') as f:
    lines = f.readlines()
    t_small= [float(line.split()[0]) for line in lines]
    z_small= [float(line.split()[1]) for line in lines]

with open('muscle_2iterations.txt', 'r') as f:
    lines = f.readlines()
    t_2iterations = [float(line.split()[0]) for line in lines]
    z_2iterations= [float(line.split()[1]) for line in lines]

with open('muscle_4iterations.txt', 'r') as f:
    lines = f.readlines()
    t_4iterations = [float(line.split()[0]) for line in lines]
    z_4iterations= [float(line.split()[1]) for line in lines]

with open('muscle_implicit.txt', 'r') as f:
    lines = f.readlines()
    t_implicit = [float(line.split()[0]) for line in lines]
    z_implicit= [float(line.split()[1]) for line in lines]

with open('muscle_case2.txt', 'r') as f:
    lines = f.readlines()
    t_case2 = [float(line.split()[0]) for line in lines]
    z_case2= [float(line.split()[1]) for line in lines]


with open('muscle_case3.txt', 'r') as f:
    lines = f.readlines()
    t_case3 = [float(line.split()[0]) for line in lines]
    z_case3= [float(line.split()[1]) for line in lines]


with open('tendon_case2.txt', 'r') as f:
    lines = f.readlines()
    t_tendon2 = [float(line.split()[0]) for line in lines]
    tendon2_begin= [float(line.split()[1]) for line in lines]
    tendon2_end= [float(line.split()[2]) for line in lines]
    tendon2_length= [float(line.split()[2])-float(line.split()[1]) for line in lines]


plt.figure(1)
plt.plot(t_low,z_low,linewidth=2,label="dt_elasticity=1.0 ms")
plt.plot(t,z,linewidth=2,label="dt_elasticity=0.1 ms")
plt.plot(t_high,z_high,linewidth=2,label="dt_elasticity=0.01 ms")

plt.title("Single muscle")
plt.xlabel("time (ms)")
plt.ylabel("muscle length (cm)")
plt.legend()
plt.show()

plt.figure(2)
plt.plot(t_high,z_high,linewidth=2,label="no tendon")
plt.plot(t_large,z_large,linewidth=2,label="large tendon")
plt.plot(t_small,z_small,linewidth=2,label="small tendon")

plt.title("Explicit coupling muscle-tendon")
plt.xlabel("time (ms)")
plt.ylabel("muscle length (cm)")
plt.legend()
plt.show()

plt.figure(3)
plt.plot(t_high,z_high,linewidth=2,label="no tendon")
plt.plot(t_small,z_small,linewidth=2,label="explicit")
plt.plot(t_2iterations,z_2iterations,linewidth=2,label="implicit (2 iterations)")
plt.plot(t_4iterations,z_4iterations,linewidth=2,label="implicit (4 iterations)")
plt.plot(t_implicit,z_implicit,linewidth=2,label="implicit")

plt.title("Implicit vs explicit coupling (small tendon)")
plt.xlabel("time (ms)")
plt.ylabel("muscle length (cm)")
plt.legend()
plt.show()

plt.figure(4)
plt.plot(t_high,z_high,linewidth=2,label="no tendon")
plt.plot(t_implicit,z_implicit,linewidth=2,label="implicit")
plt.plot(t_case2,z_case2,linewidth=2,label="implicit (case 2)")
plt.plot(t_case3,z_case3,linewidth=2,label="implicit (case 3)")


plt.title("Different acceleration schemes")
plt.xlabel("time (ms)")
plt.ylabel("muscle length (cm)")
plt.legend()
plt.show()

plt.figure(5)
plt.plot(t_tendon2,tendon2_length,linewidth=2,label="implicit (case 2)")

plt.title("Tendon length")
plt.xlabel("time (ms)")
plt.ylabel("tendon length (cm)")
plt.legend()
plt.show()

