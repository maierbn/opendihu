
# load results
import numpy as np
import glob
import scipy.integrate
import py_reader    # reader utility for opendihu *.py files

# get files in directory
filenames = [filename for filename in glob.glob("build_release/out_neumann/*.py")]

# load all files
data = py_reader.load_data(filenames)

# sort by time step
data = sorted(data, key=lambda x: x['timeStepNo'])

# diffusion constant of the problem
D = 3
size_x = 10
size_y = 10

# analytic solution

# initial values
def u_initial(x, y):
  d = np.minimum(np.sqrt((x-2)**2 + (y-size_y/2)**2), np.pi)
  return 1 + np.cos(d)

# Green function
def G(x, y, xp, yp, t):
    return np.exp(-((x-xp)**2 + (y-yp)**2) / (4*D*t)) / (4*np.pi*D*t)

def G_N(x, y, xp, yp, t):
    return G(x, y, xp, yp, t) + G(x, y, -xp, yp, t)

def analytic_solution(x, y, t):
    result, abserr = scipy.integrate.dblquad(lambda yp,xp: G_N(x,y,xp,yp,t) * u_initial(xp,yp), 0, size_x, 0, size_y)
    return result 


sample_points = [(3,3), (5,5), (5,10), (13,7), (10,10), (10,15)]
n_sample_points = len(sample_points)
numeric_list = [[] for _ in range(n_sample_points)]
analytic_list = [[] for _ in range(n_sample_points)]
t_list = []

# compute errors
for timestep_no, dataset in enumerate(data):
    
    t = dataset['currentTime']
    if t < 1e-12:
        continue

    #print(f"{timestep_no} timestep no. {dataset['timeStepNo']} t: {t}")
    dataset = dataset['data']
    numeric_solution = np.array([field_variable['components'][0]['values'] for field_variable in dataset if field_variable['name'] == 'solution'][0])

    n_solution_values = len(numeric_solution)
    nx = int(np.round(n_solution_values**(1/2)-1)/2)

    # number of nodes
    mx = 2*nx + 1
    my = 2*nx + 1
    print(mx,my)
    

    for k,sample_point in enumerate(sample_points):
        i = sample_point[0]
        j = sample_point[1]
        
        x = size_x*i/(mx-1)
        y = size_y*j/(my-1)
        analytic_solution_value = analytic_solution(x, y, t)
        
        error = (numeric_solution[j*mx + i]-analytic_solution_value)/analytic_solution_value
        numeric_list[k].append(numeric_solution[j*mx + i])
        analytic_list[k].append(analytic_solution_value)

        print(f"t: {t} {i} {j} {k} numeric: {numeric_solution[j*mx + i]:.3f}, analytic: {analytic_solution_value:.3f} ({error:.3f})")

    t_list.append(t)

import matplotlib.pyplot as plt
import matplotlib

plt.figure(figsize=(10,6))
matplotlib.rcParams.update({'font.size': 18})
for k,sample_point in enumerate(sample_points):
    p, = plt.plot(t_list, numeric_list[k], "-", label=f"{sample_point}")
    plt.plot(t_list, analytic_list[k], "--", color=p.get_color())

plt.legend()
plt.grid()
plt.xlabel("t [s]")
plt.savefig("error_neumann.png")
plt.show()
