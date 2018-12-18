// This file includes all header files that may be needed from an example
#include "utility/python_utility.h"

#include "control/dihu_context.h"
#include "control/multiple_instances.h"
#include "control/load_balancing/load_balancing.h"

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
#include "time_stepping_scheme/multidomain_solver.h"

#include "spatial_discretization/finite_element_method/05_time_stepping.h"

//#include "model_order_reduction/mor.h"
#include "model_order_reduction/time_stepping_scheme_ode_reduced_explicit.h"
//#include "model_order_reduction/time_stepping_scheme_ode_reduced_implicit.h"

#include "postprocessing/streamline_tracer.h"
#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

#include "cellml/03_cellml_adapter.h"

#include "quadrature/gauss.h"
#include "quadrature/clenshaw_curtis.h"
#include "quadrature/newton_cotes.h"
#include "quadrature/mixed.h"
