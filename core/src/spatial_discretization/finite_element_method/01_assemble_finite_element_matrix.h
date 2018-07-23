#pragma once

#include "spatial_discretization/finite_element_method/00_base.h"
#include "basis_on_mesh/mixed_basis_on_mesh.h"

namespace SpatialDiscretization
{

/**
 * Class that creates the stiffness matrix by integrating its integrand over the elements.
 * What to integrate is given by the class template Term.
 */
template<typename BasisOnMeshType, typename QuadratureType, typename Term>
class AssembleFiniteElementMatrix :
  public FiniteElementMethodBase<BasisOnMeshType, QuadratureType, Term>
{
public:
  // use constructor of base class
  using FiniteElementMethodBase<BasisOnMeshType, QuadratureType, Term>::FiniteElementMethodBase;

protected:
  //! set entries in mass matrix
  void setStiffnessMatrix();
  
  //! set entries in mass matrix
  void setMassMatrix();

};

};  // namespace

#include "spatial_discretization/finite_element_method/01_assemble_stiffness_matrix.tpp"
#include "spatial_discretization/finite_element_method/01_assemble_mass_matrix.tpp"
