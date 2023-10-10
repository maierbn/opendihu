import numpy as np
import matplotlib.pyplot as plt


with open('build_release/defaulttendon.txt', 'r') as f:
    lines = f.readlines()
    tendont = [float(line.split()[0]) for line in lines]
    tendons = [float(line.split()[1]) for line in lines]
    tendone = [float(line.split()[2]) for line in lines]

with open('build_release/defaultmuscle_left.txt', 'r') as f:
    lines = f.readlines()
    leftt = [float(line.split()[0]) for line in lines]
    lefts = [float(line.split()[1]) for line in lines]
    lefte = [float(line.split()[2]) for line in lines]

with open('build_release/defaultmuscle_right.txt', 'r') as f:
    lines = f.readlines()
    rightt = [float(line.split()[0]) for line in lines]
    rights = [float(line.split()[1]) for line in lines]
    righte = [float(line.split()[2]) for line in lines]


# paper: markersize: 0.8, font.size: 18

plt.rcParams.update({'font.size': 16})
linewidthh = 1.5

plt.figure().set_figheight(6)
plt.plot(leftt,np.array(lefte) - np.array(lefts)-15.0,linewidth=linewidthh,label="Muscle 1")
plt.plot(rightt,np.array(righte) - np.array(rights)-15.0, linewidth=linewidthh,label="Muscle 2")
plt.plot(tendont,np.array(tendone) - np.array(tendons)-5.0, linewidth=linewidthh,label="Tendon")
plt.xlabel("Time (ms)")
plt.ylabel("Change in length (cm)")
plt.title("No average - 0.0")
plt.legend()
plt.show()


