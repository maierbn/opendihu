
# load results
import numpy as np
import glob
import py_reader    # reader utility for opendihu *.py files

# Hodgkin-Huxley
# ---------------
# load all files
filenames = [filename for filename in glob.glob("build_release/out_hodgkin_huxley/*.py")]
data_opendihu = py_reader.load_data(filenames)

# sort by time step
data_opendihu = sorted(data_opendihu, key=lambda x: x['timeStepNo'])

# load result from opencor
data_opencor = np.genfromtxt("hodgkin_huxley_1952_data_istim10_opencor.csv", delimiter=",", skip_header=1)
# t, Vm, n, h, m

t_list = []
states_lists_opendihu = [[] for _ in range(4)]
states_lists_opencor = [[] for _ in range(4)]

# compute errors
for timestep_no, (dataset_opendihu, dataset_opencor) in enumerate(zip(data_opendihu, data_opencor)):
    
    t = dataset_opendihu['currentTime']

    opendihu_vm = float(py_reader.get_values(dataset_opendihu, "solution", "membrane/V"))
    opendihu_n = float(py_reader.get_values(dataset_opendihu, "solution", "potassium_channel_n_gate/n"))
    opendihu_h = float(py_reader.get_values(dataset_opendihu, "solution", "sodium_channel_h_gate/h"))
    opendihu_m = float(py_reader.get_values(dataset_opendihu, "solution", "sodium_channel_m_gate/m"))

    # store opendihu state values of current timestep 
    states_lists_opendihu[0].append(opendihu_vm)
    states_lists_opendihu[1].append(opendihu_n)
    states_lists_opendihu[2].append(opendihu_h)
    states_lists_opendihu[3].append(opendihu_m)

    # store opencor state values of current timestep
    states_lists_opencor[0].append(dataset_opencor[1])
    states_lists_opencor[1].append(dataset_opencor[2])
    states_lists_opencor[2].append(dataset_opencor[3])
    states_lists_opencor[3].append(dataset_opencor[4])

    t_list.append(t)

state_names = ["Vm", "n", "h", "m"]
for i in range(4):
    error = np.array(states_lists_opendihu[i]) - np.array(states_lists_opencor[i])
    l2_error = np.linalg.norm(error) / np.sqrt(len(states_lists_opendihu[i]))

    print(f"{state_names[i]} L2 error: {l2_error:.1e}")

import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams.update({'font.size': 16})


f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(10,6))

p, = ax1.plot(t_list, states_lists_opendihu[0], "-", label=f"Vm opendihu")
ax1.plot(t_list, states_lists_opencor[0], "--", color=p.get_color(), label="Vm opencor")

for i in range(1,4):
    p, = ax2.plot(t_list, states_lists_opendihu[i], "-", label=f"{state_names[i]} opendihu")
    ax2.plot(t_list, states_lists_opencor[i], "--", color=p.get_color(), label=f"{state_names[i]} opencor")

ax1.legend(bbox_to_anchor=(1, 1))
ax2.legend(bbox_to_anchor=(1, 1))
ax1.grid()
ax2.grid()
ax1.set_ylabel("[mV]")
ax2.set_ylabel("[-]")
ax1.set_xlabel("t [ms]")
ax2.set_xlabel("t [ms]")
plt.tight_layout()
plt.savefig("hodgkin_huxley.png")
plt.show()
