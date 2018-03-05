// This file includes all header files that may be needed from an example

#include "control/dihu_context.h"
#include "control/computation.h"

#include "basis_function/lagrange.h"
#include "basis_function/hermite.h"
#include "basis_function/mixed.h"

#include "equation/laplace.h"
#include "equation/diffusion.h"
#include "equation/reaction_diffusion.h"
#include "equation/solid_mechanics.h"

#include "mesh/structured_regular_fixed.h"
#include "mesh/structured_deformable.h"
#include "mesh/unstructured_deformable.h"

#include "operator_splitting/godunov.h"

#include "spatial_discretization/finite_element_method/finite_element_method.h"


#include "time_stepping_scheme/crank_nicholson.h"
#include "time_stepping_scheme/explicit_euler.h"

#include "model_order_reduction/pod.h"

#include "cellml/cellml_adapter.h"

#include "quadrature/gauss.h"
#include "quadrature/mixed.h"