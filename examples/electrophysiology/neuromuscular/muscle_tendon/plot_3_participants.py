import numpy as np
import matplotlib.pyplot as plt

with open('muscle_right.txt', 'r') as f:
    lines = f.readlines()
    t_right = [float(line.split()[0]) for line in lines]
    z_right = [2*14.8 + 2 - float(line.split()[1]) for line in lines]


with open('muscl_left.txt', 'r') as f:
    lines = f.readlines()
    t_left = [float(line.split()[0]) for line in lines]
    z_left= [float(line.split()[1]) for line in lines]

with open('tendon.txt', 'r') as f:
    lines = f.readlines()
    t_tendon = [float(line.split()[0]) for line in lines]
    tendon_begin= [float(line.split()[1]) for line in lines]
    tendon_end= [float(line.split()[2]) for line in lines]
    tendon_length= [float(line.split()[2])-float(line.split()[1]) for line in lines]


plt.figure(1)
plt.plot(t_left,z_left, 'o', markersize=2,label="active muscle")
plt.plot(t_right,z_right,'o', markersize=2,label="passive muscle")

plt.title("3 participants (dt = 0.01)")
plt.xlabel("time (ms)")
plt.ylabel("muscle length (cm)")
plt.legend()
plt.show()

plt.figure(2)
plt.plot(t_tendon,tendon_length,'o', markersize=2,label="tendon")

plt.title("3 participants (dt = 0.01)")
plt.xlabel("time (ms)")
plt.ylabel("tendon length (cm)")
plt.legend()
plt.show()

