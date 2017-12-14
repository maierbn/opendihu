#pragma once

#include "spatial_discretization/finite_element_method/base.h"

namespace SpatialDiscretization
{

template<typename MeshType, typename BasisFunctionType, typename IntegratorType>
class FiniteElementMethodStiffnessMatrix :
  public FiniteElementMethodBase<MeshType, BasisFunctionType, IntegratorType>
{
public:
  using FiniteElementMethodBase<MeshType, BasisFunctionType, IntegratorType>::FiniteElementMethodBase;
  
protected:
 void setStiffnessMatrix();
};

// partial specialisation for RegularMesh
template<typename IntegratorType>
class FiniteElementMethodStiffnessMatrix<Mesh::RegularFixed<1ul>, BasisFunction::Lagrange<1>, IntegratorType> :
  public FiniteElementMethodBase<Mesh::RegularFixed<1ul>, BasisFunction::Lagrange<1>, IntegratorType>
{
public:
  using FiniteElementMethodBase<Mesh::RegularFixed<1ul>, BasisFunction::Lagrange<1>, IntegratorType>::FiniteElementMethodBase;
  
protected:
 void setStiffnessMatrix();
};

template<typename IntegratorType>
class FiniteElementMethodStiffnessMatrix<Mesh::RegularFixed<2ul>, BasisFunction::Lagrange<1>, IntegratorType> :
  public FiniteElementMethodBase<Mesh::RegularFixed<2ul>, BasisFunction::Lagrange<1>, IntegratorType>
{
public:
  using FiniteElementMethodBase<Mesh::RegularFixed<2ul>, BasisFunction::Lagrange<1>, IntegratorType>::FiniteElementMethodBase;
  
protected:
 void setStiffnessMatrix();
};

template<typename IntegratorType>
class FiniteElementMethodStiffnessMatrix<Mesh::RegularFixed<3ul>, BasisFunction::Lagrange<1>, IntegratorType> :
  public FiniteElementMethodBase<Mesh::RegularFixed<3ul>, BasisFunction::Lagrange<1>, IntegratorType>
{
public:
  using FiniteElementMethodBase<Mesh::RegularFixed<3ul>, BasisFunction::Lagrange<1>, IntegratorType>::FiniteElementMethodBase;
  
protected:
 void setStiffnessMatrix();
};

template<typename BasisFunctionType, typename IntegratorType>
class FiniteElementMethodStiffnessMatrix<Mesh::Deformable<1ul>, BasisFunctionType, IntegratorType> :
  public FiniteElementMethodBase<Mesh::Deformable<1ul>, BasisFunctionType, IntegratorType>
{
public:
  using FiniteElementMethodBase<Mesh::Deformable<1ul>, BasisFunctionType, IntegratorType>::FiniteElementMethodBase;
  
protected:
 void setStiffnessMatrix();
};

template<typename BasisFunctionType, typename IntegratorType>
class FiniteElementMethodStiffnessMatrix<Mesh::Deformable<2ul>, BasisFunctionType, IntegratorType> :
  public FiniteElementMethodBase<Mesh::Deformable<2ul>, BasisFunctionType, IntegratorType>
{
public:
  using FiniteElementMethodBase<Mesh::Deformable<2ul>, BasisFunctionType, IntegratorType>::FiniteElementMethodBase;
  
protected:
 void setStiffnessMatrix();
};

template<typename BasisFunctionType, typename IntegratorType>
class FiniteElementMethodStiffnessMatrix<Mesh::Deformable<3ul>, BasisFunctionType, IntegratorType> :
  public FiniteElementMethodBase<Mesh::Deformable<3ul>, BasisFunctionType, IntegratorType>
{
public:
  using FiniteElementMethodBase<Mesh::Deformable<3ul>, BasisFunctionType, IntegratorType>::FiniteElementMethodBase;
  
protected:
 void setStiffnessMatrix();
};

 
};  // namespace