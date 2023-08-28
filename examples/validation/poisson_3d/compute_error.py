
# load results
import numpy as np
import py_reader    # reader utility for opendihu *.py files

def print_errors():

    # load all files
    data = py_reader.load_data(["build_release/out/linear_regular.py", "build_release/out/quadratic_regular.py",
                                "build_release/out/linear_structured.py", "build_release/out/quadratic_structured.py"])

    # analytic solution
    def u(x,y,z):
        return x**3*y**2*z + 4*x**2*y**2*z**3 + 2*x*y**3*z - y*z**2 + 3*x**2*y + 1

    # compute errors
    # linear scenario, regular mesh
    # -----------------------------
    linear_data = data[0]['data']
    numeric_solution = np.array([field_variable['components'][0]['values'] for field_variable in linear_data if field_variable['name'] == 'solution'][0])
    
    n_solution_values = len(numeric_solution)
    n = int((n_solution_values/24)**(1/3))

    # number of nodes
    nx = 2*n
    ny = 3*n
    nz = 4*n
    mx = nx + 1
    my = ny + 1
    mz = nz + 1
    
    analytic_solution_list = [] 
    for k in range(mz):
        for j in range(my):
            for i in range(mx):
                x = 2*i/(mx-1)
                y = 3*j/(my-1)
                z = 4*k/(mz-1)
                analytic_solution = u(x, y, z)
                analytic_solution_list.append(analytic_solution)
    analytic_solution = np.array(analytic_solution_list)

    error = numeric_solution - analytic_solution

    l2_error_linear = np.linalg.norm(error) / np.sqrt(n_solution_values-1)

    print("linear")
    print(f"{n_solution_values} solution values")
    print(f"L2 error: {l2_error_linear}")

    # quadratic scenario, regular mesh
    # ----------------------------------
    quadratic_data = data[1]['data']
    numeric_solution = np.array([field_variable['components'][0]['values'] for field_variable in quadratic_data if field_variable['name'] == 'solution'][0])
    
    n_solution_values = len(numeric_solution)

    # number of nodes
    mx = 2*nx + 1
    my = 2*ny + 1
    mz = 2*nz + 1

    analytic_solution_list = [] 
    for k in range(mz):
        for j in range(my):
            for i in range(mx):
                x = 2*i/(mx-1)
                y = 3*j/(my-1)
                z = 4*k/(mz-1)
                analytic_solution = u(x, y, z)
                analytic_solution_list.append(analytic_solution)
    analytic_solution = np.array(analytic_solution_list)

    error = numeric_solution - analytic_solution

    l2_error_quadratic = np.linalg.norm(error) / np.sqrt(n_solution_values-1)

    print("\nquadratic")
    print(f"{n_solution_values} solution values")
    print(f"L2 error: {l2_error_quadratic}")

    # linear scenario, structured mesh
    # --------------------------------
    linear_data = data[2]['data']
    numeric_solution = np.array([field_variable['components'][0]['values'] for field_variable in linear_data if field_variable['name'] == 'solution'][0])
    
    n_solution_values = len(numeric_solution)
    n = int(np.round(n_solution_values**(1/3))-1)
    print(f"n_solution_values={n_solution_values}, n={n}")

    # number of nodes
    nx = n
    ny = n
    nz = n
    mx = nx + 1
    my = ny + 1
    mz = nz + 1
    
    analytic_solution_list = [] 
    for k in range(mz):
        for j in range(my):
            for i in range(mx):
                x = 2*i/(mx-1)
                y = 3*j/(my-1)
                z = 4*k/(mz-1)
                analytic_solution = u(x, y, z)
                analytic_solution_list.append(analytic_solution)
    analytic_solution = np.array(analytic_solution_list)

    error = numeric_solution - analytic_solution

    l2_error_linear_structured = np.linalg.norm(error) / np.sqrt(n_solution_values-1)

    print("\nlinear structured")
    print(f"{n_solution_values} solution values")
    print(f"L2 error: {l2_error_linear_structured}")

    # quadratic scenario, regular mesh
    # ----------------------------------
    quadratic_data = data[3]['data']
    numeric_solution = np.array([field_variable['components'][0]['values'] for field_variable in quadratic_data if field_variable['name'] == 'solution'][0])
    
    n_solution_values = len(numeric_solution)
    n = int(np.round(n_solution_values**(1/3))/2-1)

    # number of nodes
    mx = 2*nx + 1
    my = 2*ny + 1
    mz = 2*nz + 1

    analytic_solution_list = [] 
    for k in range(mz):
        for j in range(my):
            for i in range(mx):
                x = 2*i/(mx-1)
                y = 3*j/(my-1)
                z = 4*k/(mz-1)
                analytic_solution = u(x, y, z)
                analytic_solution_list.append(analytic_solution)
    analytic_solution = np.array(analytic_solution_list)

    error = numeric_solution - analytic_solution

    l2_error_quadratic_structured = np.linalg.norm(error) / np.sqrt(n_solution_values-1)

    print("\nquadratic structured")
    print(f"{n_solution_values} solution values")
    print(f"L2 error: {l2_error_quadratic_structured}")

    return l2_error_linear, l2_error_quadratic, l2_error_linear_structured, l2_error_quadratic_structured

if __name__ == "__main__":
    print_errors()