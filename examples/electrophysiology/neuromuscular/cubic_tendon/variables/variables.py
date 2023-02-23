# scenario name for log file
scenario_name = "tendon"

# timing parameters
# -----------------
end_time = 1.0                   # [ms] end time of the simulation
dt_elasticity = 0.1               # [ms] time step width of elasticity solver

# partitioning
# ------------
# this has to match the total number of processes
n_subdomains_x = 1
n_subdomains_y = 1
n_subdomains_z = 1


# tendon geometry
# -----------------
tendon_extent = [3.0, 3.0, 2.0] # [cm, cm, cm]
tendon_offset = [0.0, 0.0, 0.0]
n_elements_tendon = [6, 6, 4] 

meshes = {}
elasticity_dirichlet_bc = {}
elasticity_neumann_bc = []


# material parameters
# --------------------
rho = 10                    # [1e-4 kg/cm^3] density of the muscle (density of water)
constant_body_force = None

tendon_material = "SaintVenantKirchoff"         #use with tendon_linear_dynamic.cpp
# tendon_material = "nonLinear"                  #use with tendon_dynamic.cpp                  


