
# load results
import numpy as np
import glob
import py_reader    # reader utility for opendihu *.py files

# Hodgkin-Huxley
# ---------------
# load all files
filenames = [filename for filename in glob.glob("build_release/out_shorten/*.py")]
data_opendihu = py_reader.load_data(filenames)

# sort by time step
data_opendihu = sorted(data_opendihu, key=lambda x: x['timeStepNo'])

print(py_reader.get_component_names(data_opendihu[0], "solution"))

opendihu_components = [
    'razumova/A_1', 'razumova/A_2', 'razumova/ATP1', 'razumova/ATP2', 'razumova/Ca_1', 'razumova/Ca_2', 'razumova/Ca_ATP1', 'razumova/Ca_ATP2', 'razumova/Ca_T_2', 'razumova/Ca_Cs1', 
    'razumova/Ca_Cs2', 'razumova/Ca_P1', 'razumova/Ca_P2', 'razumova/Ca_SR1', 'razumova/Ca_SR2', 'razumova/Ca_CaT2', 'razumova/D_0', 'razumova/D_1', 'razumova/D_2', 
    'razumova/Mg1', 'razumova/Mg2', 'razumova/Mg_ATP1', 'razumova/Mg_ATP2', 'razumova/Mg_P1', 'razumova/Mg_P2', 'razumova/P', 'razumova/P_C_SR', 'razumova/P_SR', 'razumova/x_1', 'razumova/x_2',
    'sternrios/C_0', 'sternrios/C_1', 'sternrios/C_2', 'sternrios/C_3', 'sternrios/C_4', 'sternrios/O_0', 'sternrios/O_1', 'sternrios/O_2', 'sternrios/O_3', 'sternrios/O_4', 
    'wal_environment/K_e', 'wal_environment/K_i', 'wal_environment/K_t', 'wal_environment/Na_e', 'wal_environment/Na_i', 'wal_environment/Na_t', 'sarco_DR_channel/h_K', 
    'sarco_DR_channel/n', 'sarco_Na_channel/h', 'sarco_Na_channel/m', 'sarco_Na_channel/S', 't_DR_channel/h_K_t', 't_DR_channel/n_t', 't_Na_channel/h_t', 't_Na_channel/m_t', 
    't_Na_channel/S_t', 'wal_environment/vS', 'wal_environment/vT']
opencor_components = [
    "A_1","A_2","ATP1","ATP2","Ca_1","Ca_2","Ca_ATP1","Ca_ATP2","Ca_CaT2","Ca_Cs1",
    "Ca_Cs2","Ca_P1","Ca_P2","Ca_SR1","Ca_SR2","Ca_T_2","D_0","D_1","D_2",
    "Mg1","Mg2","Mg_ATP1","Mg_ATP2",
    "Mg_P1","Mg_P2","P", "P_C_SR", "P_SR", "x_1", "x_2", 
    "C_0", "C_1", "C_2", "C_3", "C_4", 
    "O_0", "O_1", "O_2", "O_3", "O_4", 
    "K_e", "K_i", "K_t", "Na_e", "Na_i", "Na_t", "h_K", 
    "n", "h", "m", "S", "h_K_t", "n_t", "h_t", "m_t", "S_t", "vS", "vT"]

# load result from opencor
data_opencor = np.genfromtxt("shorten_ocallaghan_davidson_soboleva_2007_no_stim_data_ihh50_opencor.csv", delimiter=",", skip_header=1)
# t, Vm, n, h, m

# remove time column
data_opencor = data_opencor[:,1:]

n_states = len(opencor_components)

t_list = []
states_lists_opendihu = [[] for _ in range(n_states)]
states_lists_opencor = [[] for _ in range(n_states)]

# compute errors
for timestep_no, (dataset_opendihu, dataset_opencor) in enumerate(zip(data_opendihu, data_opencor)):
    
    t = dataset_opendihu['currentTime']

    for i in range(n_states):
        opendihu_state = float(py_reader.get_values(dataset_opendihu, "solution", opendihu_components[i]))
        opencor_state = dataset_opencor[i]

        # store opendihu state values of current timestep 
        states_lists_opendihu[i].append(opendihu_state)

        # store opencor state values of current timestep
        states_lists_opencor[i].append(opencor_state)

    t_list.append(t)

for i in range(n_states):
    error = np.array(states_lists_opendihu[i]) - np.array(states_lists_opencor[i])
    l2_error = np.linalg.norm(error) / np.sqrt(len(states_lists_opendihu[i]))

    print(f"{opencor_components[i]} L2 error: {l2_error:.1e}")

import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams.update({'font.size': 16})


f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(13,6))

p, = ax1.plot(t_list, states_lists_opendihu[n_states-1], "-", label=f"vT opendihu")
ax1.plot(t_list, states_lists_opencor[n_states-1], "--", color=p.get_color(), label="vT opencor")

for i in range(n_states):
    if i == n_states-1:
        continue
    p, = ax2.plot(t_list, states_lists_opendihu[i], "-", label=f"{opencor_components[i]}")
    ax2.plot(t_list, states_lists_opencor[i], "--", color=p.get_color())

ax1.legend(bbox_to_anchor=(1, 1))
ax2.legend(bbox_to_anchor=(1, 1), ncol=6, fontsize="9")
ax1.grid()
ax2.grid()
ax1.set_ylabel("[mV]")
ax1.set_xlabel("t [ms]")
ax2.set_xlabel("t [ms]")
plt.tight_layout()
plt.savefig("shorten.png")
plt.show()
