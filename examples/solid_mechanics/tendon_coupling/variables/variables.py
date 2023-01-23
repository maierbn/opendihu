
dt_elasticity = 0.05
end_time = 23.0

tendon_material = "nonLinear"
#tendon_material = "linear"

n_elements_single_tendon = [4,4,16]
single_tendon_extent = [2, 2, 8] #cm
single_tendon_offset = [0,0,0] #cm 

rho = 10   ## [1e-4 kg/cm^3] density of the water

force = 1e5

elasticity_dirichlet_bc = {}