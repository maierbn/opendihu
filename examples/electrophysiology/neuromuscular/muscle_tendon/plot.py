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

# with open('muscle_precice.txt', 'r') as f:
#     lines = f.readlines()
#     t_precice = [float(line.split()[0]) for line in lines]
#     z_precice= [float(line.split()[1]) for line in lines]

# with open('muscle_precice1.txt', 'r') as f:
#     lines = f.readlines()
#     t_precice1 = [float(line.split()[0]) for line in lines]
#     z_precice1= [float(line.split()[1]) for line in lines]

# with open('muscle_precice2.txt', 'r') as f:
#     lines = f.readlines()
#     t_precice2 = [float(line.split()[0]) for line in lines]
#     z_precice2= [float(line.split()[1]) for line in lines]

# with open('muscle_precice_explicit.txt', 'r') as f:
#     lines = f.readlines()
#     t_explicit = [float(line.split()[0]) for line in lines]
#     z_explicit= [float(line.split()[1]) for line in lines]


# with open('tendon_precice_explicit.txt', 'r') as f:
#     lines = f.readlines()
#     t_tendon = [float(line.split()[0]) for line in lines]
#     tendon_begin= [float(line.split()[1]) for line in lines]
#     tendon_end= [float(line.split()[2]) for line in lines]
#     tendon_length= [float(line.split()[2])-float(line.split()[1]) for line in lines]



# with open('tendon_precice1.txt', 'r') as f:
#     lines = f.readlines()
#     t_tendon1 = [float(line.split()[0]) for line in lines]
#     tendon1_begin= [float(line.split()[1]) for line in lines]
#     tendon1_end= [float(line.split()[2]) for line in lines]
#     tendon1_length= [float(line.split()[2])-float(line.split()[1]) for line in lines]

# with open('tendon_precice2.txt', 'r') as f:
#     lines = f.readlines()
#     t_tendon2 = [float(line.split()[0]) for line in lines]
#     tendon2_begin= [float(line.split()[1]) for line in lines]
#     tendon2_end= [float(line.split()[2]) for line in lines]
#     tendon2_length= [float(line.split()[2])-float(line.split()[1]) for line in lines]


plt.figure(1)
plt.plot(t_low,z_low, "o", label="dt_elasticity=1.0 ms")
plt.plot(t,z, "o", label="dt_elasticity=0.1 ms")
plt.plot(t_high,z_high, "o", label="dt_elasticity=0.01 ms")

plt.title("Single muscle")
plt.xlabel("time (ms)")
plt.ylabel("muscle length (cm)")
plt.legend()
plt.show()

# plt.figure(2)
# plt.plot(t_high,z_high, "-", label="reference")
# plt.plot(t_explicit,z_explicit, "x", label="explicit")
# plt.plot(t_precice,z_precice, "o", label="implicit v1")
# plt.plot(t_precice1,z_precice1, "o", label="implicit v2")


# plt.title("Single muscle vs coupled muscle")
# plt.xlabel("time (ms)")
# plt.ylabel("displacement (cm)")
# plt.legend()
# plt.show()

# plt.figure(3)
# #plt.plot(t_tendon,tendon_length, "o")
# plt.plot(t_tendon,tendon_end, "x", label="explicit")
# plt.plot(t_tendon1,tendon1_end, "o", label="implicit v1")
# plt.plot(t_tendon2,tendon2_end, "o", label="implicit v2")

# plt.title("Displacement of the end of the tendon")
# plt.xlabel("time (ms)")
# plt.ylabel("displacement of the free end of the tendon (cm)")
# plt.show()

# plt.figure(4)
# #plt.plot(t_tendon,tendon_length, "o")
# plt.plot(t_tendon,tendon_length, "x", label="explicit")
# plt.plot(t_tendon1,tendon1_length, "o", label="implicit v1")
# plt.plot(t_tendon2,tendon2_length, "o", label="implicit v2")

# plt.title("Length of the tendon")
# plt.xlabel("time (ms)")
# plt.ylabel("length (cm)")
# plt.show()
