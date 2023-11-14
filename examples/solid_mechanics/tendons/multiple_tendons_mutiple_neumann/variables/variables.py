# scenario name for log file
case_name = "default"
scenario_name = ""
precice_file = "../precice_config.xml"

# timing parameters
# -----------------
dt_elasticity = 1.0
end_time = 100.0

output_timestep_elasticity = 10      # [ms] timestep for elasticity output files


# partitioning
# ------------
# this has to match the total number of processes
n_subdomains_x = 1
n_subdomains_y = 1
n_subdomains_z = 1

tendon_material = "nonLinear"
#tendon_material = "SaintVenantKirchoff"         

tendon_extent = [3.0, 3.0, 5.0] # [cm, cm, cm] 2.96
n_elements_tendon = [8,8,8] 

left_tendon_offset = [0.0, 0.0, 0.0]
central_tendon_offset = [0.0, 0.0, tendon_extent[2]]
right_tendon_offset = [0.0, 0.0, 2*tendon_extent[2]]


rho = 10   ## [1e-4 kg/cm^3] density of the water

force = 1e5

elasticity_dirichlet_bc = {}
elasticity_neumann_bc = []
meshes = {}


# solvers
# -------
diffusion_solver_type = "cg"        # solver and preconditioner for the diffusion part of the Monodomain equation
diffusion_preconditioner_type = "none"      # preconditioner
diffusion_solver_maxit = 1e4
diffusion_solver_reltol = 1e-10

elasticity_solver_type = "preonly"
elasticity_preconditioner_type = "lu"
snes_max_iterations = 10                  # maximum number of iterations in the nonlinear solver
snes_rebuild_jacobian_frequency = 2       # how often the jacobian should be recomputed, -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time the Jacobian is built etc. -2 means rebuild at next chance but then never again 
snes_relative_tolerance = 1e-5      # relative tolerance of the nonlinear solver
snes_absolute_tolerance = 1e-5      # absolute tolerance of the nonlinear solver
linear_relative_tolerance = 1e-5           # relative tolerance of the residual of the linear solver
linear_absolute_tolerance = 1e-10          # absolute tolerance of the residual of the linear solver

optimization_type = "vc"            # the optimization_type used in the cellml adapter, "vc" uses explicit vectorization
approximate_exponential_function = False   # if the exponential function should be approximated by a Taylor series with only 11 FLOPS
dynamic = True                      # if the dynamic hyperelasticity solver should be used

def get_from_obj(data, path):
    for elem in path:
        if type(elem) == str:
            data = data[elem]
        elif type(elem) == int:
            data = data[elem]
        elif type(elem) == tuple:
            # search for key == value with (key, value) = elem
            key, value = elem
            data = next(filter(lambda e: e[key] == value, data))
        else:
            raise KeyError(f"Unknown type of '{elem}': '{type(elem)}'. Path: '{'.'.join(path)}'")
    return data

# available slots for data fields

# geometry
# u
# v
# PK2-Stress (Voigt)
# active PK2-Stress (Voigt)
# fiberDirection
# t (current traction)
# T (material traction)
# F
# Fdot
# P (PK1 stress)
# σ (Cauchy stress)
# J

def muscle_left_write_to_file(data):
    t = get_from_obj(data, [0, 'currentTime'])
    z_data = get_from_obj(data, [0, 'data', ('name','geometry'), 'components', 2, 'values'])

    [mx, my, mz] = get_from_obj(data, [0, 'nElementsLocal'])
    nx = 2*mx + 1
    ny = 2*my + 1
    nz = 2*mz + 1
    # compute average z-value of end of muscle
    z_value_begin = 0.0
    z_value_end = 0.0

    for j in range(ny):
        for i in range(nx):
            z_value_begin += z_data[j*nx + i]
            z_value_end += z_data[(nz-1)*nx*ny + j*nx + i]


    z_value_begin /= ny*nx
    z_value_end /= ny*nx


    f = open(case_name + "muscle_left.txt", "a")
    f.write("{:6.2f} {:+2.8f} {:+2.8f}\n".format(t,z_value_begin, z_value_end))
    f.close()

def muscle_right_write_to_file(data):
    t = get_from_obj(data, [0, 'currentTime'])
    z_data = get_from_obj(data, [0, 'data', ('name','geometry'), 'components', 2, 'values'])

    [mx, my, mz] = get_from_obj(data, [0, 'nElementsLocal'])
    nx = 2*mx + 1
    ny = 2*my + 1
    nz = 2*mz + 1
    # compute average z-value of end of muscle
    z_value_begin = 0.0
    z_value_end = 0.0

    for j in range(ny):
        for i in range(nx):
            z_value_begin += z_data[j*nx + i]
            z_value_end += z_data[(nz-1)*nx*ny + j*nx + i]

    z_value_begin /= ny*nx
    z_value_end /= ny*nx

    f = open(case_name + "muscle_right.txt", "a")
    f.write("{:6.2f} {:+2.8f} {:+2.8f}\n".format(t,z_value_begin, z_value_end))
    f.close()

def tendon_write_to_file(data):
    t = get_from_obj(data, [0, 'currentTime'])
    z_data = get_from_obj(data, [0, 'data', ('name','geometry'), 'components', 2, 'values'])

    [mx, my, mz] = get_from_obj(data, [0, 'nElementsLocal'])
    nx = 2*mx + 1
    ny = 2*my + 1
    nz = 2*mz + 1
    # compute average z-value of end of muscle
    z_value_begin = 0.0
    z_value_end = 0.0

    for j in range(ny):
        for i in range(nx):
            z_value_begin += z_data[j*nx + i]
            z_value_end += z_data[(nz-1)*nx*ny + j*nx + i]

    z_value_begin /= ny*nx
    z_value_end /= ny*nx


    f = open(case_name + "tendon.txt", "a")
    f.write("{:6.2f} {:+2.8f} {:+2.8f}\n".format(t,z_value_begin, z_value_end))
    f.close()

