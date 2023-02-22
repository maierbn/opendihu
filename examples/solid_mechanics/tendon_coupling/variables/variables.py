dt_elasticity = 0.01
end_time = 5.0

tendon_material = "nonLinear"
#tendon_material = "linear"

n_elements_single_tendon = [4,4,16]
single_tendon_extent = [2, 2, 8] #cm
single_tendon_offset = [0,0,0] #cm 

rho = 10   ## [1e-4 kg/cm^3] density of the water

force = 1e5

elasticity_dirichlet_bc = {}

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

def tendon_postprocess(data):
    t = get_from_obj(data, [0, 'currentTime'])
    z_data = get_from_obj(data, [0, 'data', ('name','geometry'), 'components', 2, 'values'])
    [mx, my, mz] = get_from_obj(data, [0, 'nElementsLocal'])
    nx = 2*mx + 1
    ny = 2*my + 1
    nz = 2*mz + 1
    # compute average z-value of end of muscle
    z_value = 0
    for j in range(ny):
        for i in range(nx):
            z_value += z_data[(nz-1)*nx*ny + j*nx + i]
    z_value /= ny*nx

    print("Tendon: t= {:6.2f}, z-position of the right side = {:+2.4f}".format(t, z_value))

def write_to_file(data):
    t = get_from_obj(data, [0, 'currentTime'])
    z_data = get_from_obj(data, [0, 'data', ('name','geometry'), 'components', 2, 'values'])
    [mx, my, mz] = get_from_obj(data, [0, 'nElementsLocal'])
    nx = 2*mx + 1
    ny = 2*my + 1
    nz = 2*mz + 1
    # compute average z-value of end of muscle
    z_value = 0
    for j in range(ny):
        for i in range(nx):
            z_value += z_data[(nz-1)*nx*ny + j*nx + i]
    z_value /= ny*nx

    f = open("displacement6.txt", "a")
    f.write("{:6.2f} {:+2.8f}\n".format(t, z_value))
    f.close()


