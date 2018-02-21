#pragma once

#include "spatial_discretization/finite_element_method/02_stiffness_matrix.h"

namespace SpatialDiscretization
{

/**
 * Class that sets the right hand side vector by integrating the integrand over the elements.
 * What to integrate is given by the class template Term.
 */
template<typename BasisOnMeshType, typename QuadratureType, typename Term, typename Dummy= Term>
class AssembleRightHandSide :
  public FiniteElementMethodStiffnessMatrix<BasisOnMeshType, QuadratureType, Term>
{
public:
  // use constructor of base class
  using FiniteElementMethodStiffnessMatrix<BasisOnMeshType, QuadratureType, Term>::FiniteElementMethodStiffnessMatrix;

protected:
  //! Transform values in rhs vector into FEM discretized values by multiplying them with the integrate basis functions
  void transferRhsToWeakForm();
  
  //! set the matrix that transforms a vector of rhs value into a vector that contains the right hand side in discretized form for FEM
  void setMassMatrix();
};
 
/**
 * Partial specialization for solid mechanics
 */
template<typename MixedBasisOnMeshType, typename MixedQuadratureType, typename Term>
class AssembleRightHandSide<MixedBasisOnMeshType, MixedQuadratureType, Term, Equation::isSolidMechanics<Term>> :
  public FiniteElementMethodStiffnessMatrix<MixedBasisOnMeshType, MixedQuadratureType, Equation::Static::SolidMechanics>
{
public:
  // use constructor of base class
  using FiniteElementMethodStiffnessMatrix<MixedBasisOnMeshType, MixedQuadratureType, Equation::Static::SolidMechanics>::FiniteElementMethodStiffnessMatrix;

protected:
  //! Transform values in rhs vector into FEM discretized values by multiplying them with the integrate basis functions
  void transferRhsToWeakForm(){}
  
  //! set the matrix that transforms a vector of rhs value into a vector that contains the right hand side in discretized form for FEM
  void setMassMatrix(){}
};

};  // namespace

#include "spatial_discretization/finite_element_method/03_assemble_rhs.tpp"
