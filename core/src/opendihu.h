// This file includes all header files that may be needed from an example

#include "control/dihu_context.h"
#include "control/computation.h"

#include "basis_function/lagrange.h"

#include "equation/laplace.h"
#include "equation/diffusion.h"
#include "equation/reaction_diffusion.h"

#include "mesh/deformable.h"
#include "mesh/nonrectilinear_fixed.h"
#include "mesh/rectilinear_fixed.h"
#include "mesh/regular_fixed.h"

#include "operator_splitting/godunov.h"

#include "spatial_discretization/finite_element_method.h"

#include "time_stepping_scheme/crank_nicholson.h"
#include "time_stepping_scheme/explicit_euler.h"

#include "cellml/cellml_adapter.h"