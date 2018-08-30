#include "spatial_discretization/finite_element_method/solid_mechanics/02_stiffness_matrix_incompressible.h"

#include <Python.h>  // has to be the first included header
#include <memory>
#include <petscksp.h>
#include <petscsys.h>

#include "semt/Semt.h"
#include "semt/Shortcuts.h"

namespace SpatialDiscretization
{

// set stiffness matrix for a u-p mixed formulation in which the pressure is condensed out
template<typename HighOrderBasisOnMeshType, int completePolynomialOrder, typename MixedQuadratureType, typename Term>
void FiniteElementMethodMatrix<
  MixedBasisOnMeshTemplate<HighOrderBasisOnMeshType, completePolynomialOrder>,
  MixedQuadratureType,
  Term,
  Mesh::isDeformable<typename HighOrderBasisOnMeshType::Mesh>,
  Equation::isIncompressible<Term>
>::
setStiffnessMatrix()
{
}

};    // namespace