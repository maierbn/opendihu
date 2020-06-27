# This file contains internal values of create_partitioned_meshes_for_settings.py 
# that need to be persistent for the helper function to work

debug_output = False                # verbose output in this python script, for debugging the domain decomposition

dt_elasticity = 1e-3

# further internal variables
n_fibers_total = None
n_subdomains_xy = None
own_subdomain_coordinate_x = None
own_subdomain_coordinate_y = None
own_subdomain_coordinate_z = None
n_fibers_x = None
n_fibers_y = None
n_points_whole_fiber = None
n_points_3D_mesh_global_x = None
n_points_3D_mesh_global_y = None
n_points_3D_mesh_global_z = None
nodal_stimulation_current = None
fiber_file_handle = None
fibers = None
fiber_distribution = None
firing_times = None
n_fibers_per_subdomain_x = None
n_fibers_per_subdomain_y = None
n_points_per_subdomain_z = None
z_point_index_start = None
z_point_index_end = None
n_elements_3D_mesh = None
meshes = None
fibers_on_own_rank = None
n_fiber_nodes_on_subdomain = None
fiber_start_node_no = None
cellml_file = ""
n_subdomains_x = 1
n_subdomains_y = 1
n_subdomains_z = 1
