import numpy as np
import matplotlib.pyplot as plt

with open('data/muscle_case4.txt', 'r') as f:
    lines = f.readlines()
    t4 = [float(line.split()[0]) for line in lines]
    z4 = [float(line.split()[1]) for line in lines]


with open('data/muscle_case5.txt', 'r') as f:
    lines = f.readlines()
    t5 = [float(line.split()[0]) for line in lines]
    z5= [float(line.split()[1]) for line in lines]

with open('data/tendon_case4.txt', 'r') as f:
    lines = f.readlines()
    t_tendon4 = [float(line.split()[0]) for line in lines]
    tendon_length4= [float(line.split()[2])-float(line.split()[1]) for line in lines]

with open('data/tendon_case5.txt', 'r') as f:
    lines = f.readlines()
    t_tendon5 = [float(line.split()[0]) for line in lines]
    tendon_length5= [float(line.split()[2])-float(line.split()[1]) for line in lines]


plt.figure(1)
plt.plot(t4,z4, 'o', markersize=2,label="direct transfer")
plt.plot(t5,z5,'o', markersize=2,label="reversed transfer")

plt.title("direct vs reversed transfer of data")
plt.xlabel("time (ms)")
plt.ylabel("muscle length (cm)")
plt.legend()
plt.show()

plt.figure(2)
plt.plot(t_tendon4,tendon_length4,'o', markersize=2,label="direct transfer")
plt.plot(t_tendon5,tendon_length5,'o', markersize=2,label="reversed transfer")

plt.title("direct vs reversed transfer of data")
plt.xlabel("time (ms)")
plt.ylabel("tendon length (cm)")
plt.legend()
plt.show()

