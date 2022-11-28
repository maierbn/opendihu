// This file includes all header files that may be needed from an example
#include "utility/python_utility.h"

#include "control/coupling/coupling.h"
#include "control/coupling/multiple_coupling.h"
#include "control/dihu_context.h"
#include "control/multiple_instances.h"
#include "control/load_balancing/load_balancing.h"
#include "control/map_dofs/map_dofs.h"
#include "control/precice/old/partitioned_fibers.h"
#include "control/precice/old/muscle_contraction.h"
#include "control/precice/old/contraction_dirichlet_boundary_conditions.h"
#include "control/precice/old/contraction_neumann_boundary_conditions.h"
#include "control/precice/surface_coupling/precice_adapter.h"
#include "control/precice/volume_coupling/precice_adapter_volume_coupling.h"

#include "basis_function/lagrange.h"
#include "basis_function/hermite.h"
#include "basis_function/mixed.h"

#include "equation/equations.h"

#include "mesh/structured_regular_fixed.h"
#include "mesh/structured_deformable.h"
#include "mesh/unstructured_deformable.h"

#include "operator_splitting/godunov.h"
#include "operator_splitting/strang.h"

#include "spatial_discretization/finite_element_method/finite_element_method.h"

#include "interfaces/discretizable_in_time.h"

#include "time_stepping_scheme/crank_nicolson.h"
#include "time_stepping_scheme/explicit_euler.h"
#include "time_stepping_scheme/implicit_euler.h"
#include "time_stepping_scheme/heun.h"
#include "time_stepping_scheme/repeated_call.h"
#include "time_stepping_scheme/repeated_call_static.h"
#include "specialized_solver/multidomain_solver/multidomain_solver.h"
#include "specialized_solver/multidomain_solver/multidomain_with_fat_solver.h"
#include "specialized_solver/static_bidomain_solver.h"
#include "specialized_solver/solid_mechanics/quasi_static/quasi_static_linear_elasticity_solver.h"
#include "specialized_solver/solid_mechanics/quasi_static/quasi_static_nonlinear_elasticity_solver_chaste.h"
#include "specialized_solver/solid_mechanics/quasi_static/febio/quasi_static_nonlinear_elasticity_solver_febio.h"
#include "specialized_solver/solid_mechanics/quasi_static/quasi_static_linear_elasticity_solver.h"
#include "specialized_solver/solid_mechanics/hyperelasticity/02_hyperelasticity_solver.h"
#include "specialized_solver/solid_mechanics/dynamic_hyperelasticity/dynamic_hyperelasticity_solver.h"
#include "specialized_solver/solid_mechanics/quasistatic_hyperelasticity/quasistatic_hyperelasticity_solver.h"
#include "specialized_solver/fast_monodomain_solver/fast_monodomain_solver.h"
#include "specialized_solver/prescribed_values.h"
#include "specialized_solver/my_new_solver/my_new_static_solver.h"
#include "specialized_solver/my_new_solver/my_new_timestepping_solver.h"
#include "specialized_solver/dummy.h"
#include "specialized_solver/muscle_contraction_solver.h"
#include "time_stepping_scheme/heun_adaptive.h"

#include "spatial_discretization/finite_element_method/05_time_stepping.h"

//#include "model_order_reduction/mor.h"
//#include "model_order_reduction/time_stepping_reduced_explicit.h"
//#include "model_order_reduction/time_stepping_reduced_implicit.h"
#include "model_order_reduction/explicit_euler_reduced.h"
#include "model_order_reduction/implicit_euler_reduced.h"

#include "postprocessing/streamline_tracer.h"
#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"
#include "postprocessing/scale_fibers_in_file.h"

#include "cellml/03_cellml_adapter.h"

#include "output_writer/output_surface/output_surface.h"

#include "quadrature/gauss.h"
#include "quadrature/clenshaw_curtis.h"
#include "quadrature/newton_cotes.h"
#include "quadrature/mixed.h"
