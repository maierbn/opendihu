#pragma once

#include "spatial_discretization/finite_element_method/00_base.h"
#include "basis_on_mesh/mixed_basis_on_mesh.h"

namespace SpatialDiscretization
{

/**
 * Class that creates the stiffness matrix by integrating its integrand over the elements.
 * What to integrate is given by the class template Term.
 */
template<typename BasisOnMeshType, typename IntegratorType, typename Term>
class AssembleStiffnessMatrix :
  public FiniteElementMethodBase<BasisOnMeshType, IntegratorType>
{
public:
  // use constructor of base class
  using FiniteElementMethodBase<BasisOnMeshType, IntegratorType>::FiniteElementMethodBase;
  
protected:
 void setStiffnessMatrix();
  
};

 
};  // namespace

#include "spatial_discretization/finite_element_method/01_assemble_stiffness_matrix.tpp"
