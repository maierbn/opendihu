import os
import subprocess
import numpy as np
import compute_error
import matplotlib.pyplot as plt
import matplotlib

errors_linear = []
errors_quadratic = []
errors_hermite = []
n_list = []

for n in np.logspace(1, 3, 10):
    n = int(n)
    n_list.append(n)
    
    subprocess.call(f"cd build_release; ./linear ../settings_poisson_1d.py {n}", shell=True)
    subprocess.call(f"cd build_release; ./quadratic ../settings_poisson_1d.py {n}", shell=True)
    subprocess.call(f"cd build_release; ./hermite ../settings_poisson_1d.py {n}", shell=True)
    l2_error_linear, l2_error_quadratic, l2_error_hermite = compute_error.print_errors()

    errors_linear.append(l2_error_linear)
    errors_quadratic.append(l2_error_quadratic)
    errors_hermite.append(l2_error_hermite)

print(f"n_list: {n_list}")
print(f"errors_linear: {errors_linear}")
print(f"errors_quadratic: {errors_quadratic}")
print(f"errors_hermite: {errors_hermite}")

print(f"experimental order of convergence (linear): {[(np.log(errors_linear[i])-np.log(errors_linear[i-1]))/(np.log(n_list[i])-np.log(n_list[i-1])) for i in range(1,len(n_list))]}")
print(f"experimental order of convergence (quadratic): {[(np.log(errors_quadratic[i])-np.log(errors_quadratic[i-1]))/(np.log(n_list[i])-np.log(n_list[i-1])) for i in range(1,len(n_list))]}")
print(f"experimental order of convergence (hermite): {[(np.log(errors_hermite[i])-np.log(errors_hermite[i-1]))/(np.log(n_list[i])-np.log(n_list[i-1])) for i in range(1,len(n_list))]}")

matplotlib.rcParams.update({'font.size': 18})
plt.figure(figsize=(10,6))
plt.plot(n_list, errors_hermite, "o-", label="hermite")
plt.plot(n_list, errors_linear, "o-", label="linear")
plt.plot(n_list, errors_quadratic, "o-", label="quadratic")

plt.xscale('log')
plt.yscale('log')
plt.grid(which='both')
plt.xlabel("number of elements")
plt.ylabel("L2 error")
plt.legend()
plt.tight_layout()

plt.savefig("convergence_study.pdf")
plt.show()
