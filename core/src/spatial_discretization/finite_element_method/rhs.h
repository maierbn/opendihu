#pragma once

#include "spatial_discretization/finite_element_method/stiffness_matrix.h"

namespace SpatialDiscretization
{
 
template<typename MeshType, typename BasisFunctionType, typename IntegratorType>
class FiniteElementMethodBaseRhs :
  public FiniteElementMethodStiffnessMatrix<MeshType, BasisFunctionType, IntegratorType>
{
};

template<typename IntegratorType>
class FiniteElementMethodBaseRhs<Mesh::RegularFixed<1ul>, BasisFunction::Lagrange<1>, IntegratorType> :
  public FiniteElementMethodStiffnessMatrix<Mesh::RegularFixed<1ul>, BasisFunction::Lagrange<1>, IntegratorType>
{
public:
  using FiniteElementMethodStiffnessMatrix<Mesh::RegularFixed<1ul>, BasisFunction::Lagrange<1>, IntegratorType>::FiniteElementMethodStiffnessMatrix;
protected:
  //! Transform values in rhs vector into FEM discretized values by multiplying them with the integrate basis functions
  void transferRhsToWeakForm();
};

template<typename IntegratorType>
class FiniteElementMethodBaseRhs<Mesh::RegularFixed<2ul>, BasisFunction::Lagrange<1>, IntegratorType> :
  public FiniteElementMethodStiffnessMatrix<Mesh::RegularFixed<2ul>, BasisFunction::Lagrange<1>, IntegratorType>
{
public:
  using FiniteElementMethodStiffnessMatrix<Mesh::RegularFixed<2ul>, BasisFunction::Lagrange<1>, IntegratorType>::FiniteElementMethodStiffnessMatrix;
protected:
  //! Transform values in rhs vector into FEM discretized values by multiplying them with the integrate basis functions
  void transferRhsToWeakForm();
};

template<typename IntegratorType>
class FiniteElementMethodBaseRhs<Mesh::RegularFixed<3ul>, BasisFunction::Lagrange<1>, IntegratorType> :
  public FiniteElementMethodStiffnessMatrix<Mesh::RegularFixed<3ul>, BasisFunction::Lagrange<1>, IntegratorType>
{
public:
  using FiniteElementMethodStiffnessMatrix<Mesh::RegularFixed<3ul>, BasisFunction::Lagrange<1>, IntegratorType>::FiniteElementMethodStiffnessMatrix;
protected:
  //! Transform values in rhs vector into FEM discretized values by multiplying them with the integrate basis functions
  void transferRhsToWeakForm();
};

template<typename BasisFunctionType, typename IntegratorType>
class FiniteElementMethodBaseRhs<Mesh::Deformable<1ul>, BasisFunctionType, IntegratorType> :
  public FiniteElementMethodStiffnessMatrix<Mesh::Deformable<1ul>, BasisFunctionType, IntegratorType>
{
public:
  using FiniteElementMethodStiffnessMatrix<Mesh::Deformable<1ul>, BasisFunctionType, IntegratorType>::FiniteElementMethodStiffnessMatrix;
protected:
  //! Transform values in rhs vector into FEM discretized values by multiplying them with the integrate basis functions
  void transferRhsToWeakForm();
};

template<typename BasisFunctionType, typename IntegratorType>
class FiniteElementMethodBaseRhs<Mesh::Deformable<2ul>, BasisFunctionType, IntegratorType> :
  public FiniteElementMethodStiffnessMatrix<Mesh::Deformable<2ul>, BasisFunctionType, IntegratorType>
{
public:
  using FiniteElementMethodStiffnessMatrix<Mesh::Deformable<2ul>, BasisFunctionType, IntegratorType>::FiniteElementMethodStiffnessMatrix;
protected:
  //! Transform values in rhs vector into FEM discretized values by multiplying them with the integrate basis functions
  void transferRhsToWeakForm();
};

template<typename BasisFunctionType, typename IntegratorType>
class FiniteElementMethodBaseRhs<Mesh::Deformable<3ul>, BasisFunctionType, IntegratorType> :
  public FiniteElementMethodStiffnessMatrix<Mesh::Deformable<3ul>, BasisFunctionType, IntegratorType>
{
public:
  using FiniteElementMethodStiffnessMatrix<Mesh::Deformable<3ul>, BasisFunctionType, IntegratorType>::FiniteElementMethodStiffnessMatrix;
protected:
  //! Transform values in rhs vector into FEM discretized values by multiplying them with the integrate basis functions
  void transferRhsToWeakForm();
};


// base class implementing right hand side, that can be set by user for poisson equation
template<typename MeshType, typename BasisFunctionType, typename IntegratorType>
class FiniteElementMethodRhs :
  public FiniteElementMethodBaseRhs<MeshType, BasisFunctionType, IntegratorType>
  //public FiniteElementMethodStiffnessMatrix<MeshType, BasisFunctionType, IntegratorType>
{
public:
  using FiniteElementMethodBaseRhs<MeshType, BasisFunctionType, IntegratorType>::FiniteElementMethodBaseRhs;
 
  friend class StiffnessMatrixTester;    ///< a class used for testing 
protected:
 
  //! read in rhs values from config and transfer to weak form
  void setRightHandSide();
};
 
};  // namespace