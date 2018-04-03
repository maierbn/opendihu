#include "spatial_discretization/finite_element_method/solid_mechanics/02_stiffness_matrix_incompressible.h"

#include <Python.h>  // has to be the first included header
#include <memory>
#include <petscksp.h>
#include <petscsys.h>

#include "semt/Semt.h"
#include "semt/Shortcuts.h"

namespace SpatialDiscretization
{
  
// general implementation for solid mechanics
template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  Mesh::isDeformable<typename BasisOnMeshType::Mesh>,
  Equation::isIncompressible<Term>
>::
setStiffnessMatrix()
{
}

};    // namespace