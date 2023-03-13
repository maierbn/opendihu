import numpy as np
import matplotlib.pyplot as plt

with open('case14/muscle.txt', 'r') as f:
    lines = f.readlines()
    timewindow_muscle = [float(line.split()[0]) for line in lines]
    muscle_end_z = [float(line.split()[1]) for line in lines]
    muscle_end_traction = [float(line.split()[2]) for line in lines]

with open('case14/tendon.txt', 'r') as f:
    lines = f.readlines()
    timewindow_tendon = [float(line.split()[0]) for line in lines]
    tendon_start_z = [float(line.split()[1]) for line in lines]
    tendon_start_traction = [float(line.split()[3]) for line in lines]


plt.figure(1)
plt.plot(timewindow_muscle,muscle_end_z,linewidth=1.5, label="muscle end")
plt.plot(timewindow_tendon,tendon_start_z, linewidth=1.5, label="tendon start")
plt.title("Case 4")
plt.xlabel("time (ms)")
plt.ylabel("z (cm)")
plt.legend()
plt.show()

plt.figure(2)
plt.plot(timewindow_muscle,muscle_end_traction, linewidth=2, label="muscle end")
plt.plot(timewindow_tendon,tendon_start_traction, linewidth=2, label="tendon start")
plt.title("Case 4")
plt.xlabel("time (ms)")
plt.ylabel("traction")
plt.legend()
plt.show()

plt.figure(3)
plt.plot(timewindow_muscle,muscle_end_traction, linewidth=2, label="muscle end")
plt.plot(timewindow_tendon,np.array(tendon_start_traction)*(-0.1), linewidth=2, label="tendon start scaled by -0.1")
plt.title("Case 4")
plt.xlabel("time (ms)")
plt.ylabel("traction")
plt.legend()
plt.show()