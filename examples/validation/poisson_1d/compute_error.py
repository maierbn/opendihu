
# load results
import numpy as np
import py_reader    # reader utility for opendihu *.py files

def print_errors():

    # load all files
    data = py_reader.load_data(["build_release/out/linear.py", "build_release/out/quadratic.py", "build_release/out/hermite.py"])

    # analytic solution
    def u(x):
        return -1/12*x**4 + 1/2*x**2 + 13/12*x + 1

    # compute errors
    # linear scenario
    # ---------------
    linear_data = data[0]['data']
    numeric_solution = np.array([field_variable['components'][0]['values'] for field_variable in linear_data if field_variable['name'] == 'solution'][0])
    n_nodes = len(numeric_solution)

    analytic_solution = u(np.linspace(0, 3, n_nodes))
    error = numeric_solution - analytic_solution

    l2_error_linear = np.linalg.norm(error) / np.sqrt(n_nodes-1)

    print("linear")
    print(f"analytic solution: {analytic_solution}")
    print(f"numeric solution: {numeric_solution}")
    print(f"error: {error}")
    print(f"L2 error: {l2_error_linear}")

    # quadratic scenario
    # ---------------
    quadratic_data = data[1]['data']
    numeric_solution = np.array([field_variable['components'][0]['values'] for field_variable in quadratic_data if field_variable['name'] == 'solution'][0])
    n_nodes = len(numeric_solution)

    analytic_solution = u(np.linspace(0, 3, n_nodes))
    error = numeric_solution - analytic_solution

    l2_error_quadratic = np.linalg.norm(error) / np.sqrt(n_nodes-1)

    print("\nquadratic")
    print(f"analytic solution: {analytic_solution}")
    print(f"numeric solution: {numeric_solution}")
    print(f"error: {error}")
    print(f"L2 error: {l2_error_quadratic}")

    # hermite scenario
    # ---------------
    hermite_data = data[2]['data']
    numeric_solution = [field_variable['components'][0]['values'] for field_variable in hermite_data if field_variable['name'] == 'solution'][0]
    numeric_solution = np.array(numeric_solution[::2])
    n_nodes = len(numeric_solution)

    analytic_solution = u(np.linspace(0, 3, n_nodes))
    error = numeric_solution - analytic_solution

    l2_error_hermite = np.linalg.norm(error) / np.sqrt(n_nodes-1)

    print("\nHermite")
    print(f"analytic solution: {analytic_solution}")
    print(f"numeric solution: {numeric_solution}")
    print(f"error: {error}")
    print(f"L2 error: {l2_error_hermite}")

    return l2_error_linear, l2_error_quadratic, l2_error_hermite

if __name__ == "__main__":
    print_errors()