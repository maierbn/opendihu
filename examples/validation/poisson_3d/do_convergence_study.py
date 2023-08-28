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

if False:
    for n in np.logspace(0.2, 0.4, 2):
        n = int(n)
        n_list.append(n)
        print(f"\nn={n}")
        
        subprocess.call(f"cd build_release; ./linear_regular ../settings_poisson_3d.py {n}", shell=True)
        subprocess.call(f"cd build_release; ./quadratic_regular ../settings_poisson_3d.py {n}", shell=True)
        l2_error_linear, l2_error_quadratic, l2_error_linear_structured, l2_error_quadratic_structured = compute_error.print_errors()

        errors_linear.append(l2_error_linear)
        errors_quadratic.append(l2_error_quadratic)

    print(f"n_list: {n_list}")
    print(f"errors_linear: {errors_linear}")
    print(f"errors_quadratic: {errors_quadratic}")

    matplotlib.rcParams.update({'font.size': 18})
    plt.figure(figsize=(10,6))
    plt.plot(n_list, errors_linear, "o-", label="linear")
    plt.plot(n_list, errors_quadratic, "o-", label="quadratic")

    plt.xscale('log')
    plt.yscale('log')
    plt.grid(which='both')
    plt.xlabel("number of elements")
    plt.ylabel("L2 error")
    plt.legend()
    plt.tight_layout()

    plt.savefig("poisson_3d_convergence_study_regular.pdf")
    plt.show()

# study with non-cubic elements (structured mesh)
for n in np.logspace(0.4, 1.3, 5):
    n = int(n)
    n_list.append(n)
    
    subprocess.call(f"cd build_release; ./linear_structured ../settings_poisson_3d.py {n}", shell=True)
    subprocess.call(f"cd build_release; ./quadratic_structured ../settings_poisson_3d.py {n}", shell=True)
    l2_error_linear, l2_error_quadratic, l2_error_linear_structured, l2_error_quadratic_structured = compute_error.print_errors()

    errors_linear.append(l2_error_linear_structured)
    errors_quadratic.append(l2_error_quadratic_structured)

print(f"n_list: {n_list}")
print(f"errors_linear: {errors_linear}")
print(f"errors_quadratic: {errors_quadratic}")

matplotlib.rcParams.update({'font.size': 18})
plt.figure(figsize=(10,6))
plt.plot(n_list, errors_linear, "o-", label="linear")
plt.plot(n_list, errors_quadratic, "o-", label="quadratic")

plt.xscale('log')
plt.yscale('log')
plt.grid(which='both')
plt.xlabel("number of elements")
plt.ylabel("L2 error")
plt.legend()
plt.tight_layout()

plt.savefig("poisson_3d_convergence_study_structured.pdf")
plt.show()
