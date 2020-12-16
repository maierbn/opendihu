This directory contains python scripts that are used to work with the muscle geometry and generate meshes.
There are scripts that can be executed from the command line:
* `create_spline_surface.py` - This script needs a STL file of a tubular muscle surface and then creates a NURBS surface approximation.
* `slice_stl.py` - From a tubular surface as STL file, extract a part between two z planes, output as enclosed STL surface.
* `plot_quadrangulation_schemes.py` - This simply plots the different triangulation schemes and can be seen as a reference implementation how to construct them.

The other scripts are modules that can be loaded from other python files:
* `stl_create_rings.py` - This contains functionality to "slice" a geometry, i.e., to extract rings/loops of a tubular surface.
  * imports `stl_create_mesh.py`
  * imports `spline_surface.py`
  
* `stl_create_mesh.py` - This contains the Python parts of the algorithms to create a 2D quadrangulation of a slice (algorithm with harmonic maps) and create a 3D hex mesh out of these slices. The substeps are implemented in the imported scripts.
  * imports `spline_surface.py`
  * imports `create_planar_mesh_helper.py`
  * imports `create_slice_tringulation.py`
  * imports `solve_laplace_problem.py`
  * imports `map_quadrangulation_to_world_space.py`
  * imports `fix_and_smooth_mesh.py`

These modules are called, e.g., by the C++ code in `examples/fiber_tracing/parallel_fiber_estimation`.
