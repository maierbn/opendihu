#pragma once

#include "spatial_discretization/finite_element_method/stiffness_matrix.h"

#include "mesh/mesh.h"
#include "time_stepping_scheme/discretizable_in_time.h"

namespace SpatialDiscretization
{
 
template<typename MeshType, typename BasisFunctionType, typename IntegratorType>
class FiniteElementMethodBaseTimeStepping : 
  public FiniteElementMethodStiffnessMatrix<MeshType, BasisFunctionType, IntegratorType>
{
public:
  using FiniteElementMethodStiffnessMatrix<MeshType, BasisFunctionType, IntegratorType>::FiniteElementMethodStiffnessMatrix;
protected:
  //! create the discretization matrix which is the mapping between strong formulated and weak formulated rhs vector
  void createRhsDiscretizationMatrix();
};

template<typename IntegratorType>
class FiniteElementMethodBaseTimeStepping<Mesh::RegularFixed<1ul>, BasisFunction::Lagrange<1>, IntegratorType> : 
  public FiniteElementMethodStiffnessMatrix<Mesh::RegularFixed<1ul>, BasisFunction::Lagrange<1>, IntegratorType>
{
public:
  using FiniteElementMethodStiffnessMatrix<Mesh::RegularFixed<1ul>, BasisFunction::Lagrange<1>, IntegratorType>::FiniteElementMethodStiffnessMatrix;

protected:
  //! create the discretization matrix which is the mapping between strong formulated and weak formulated rhs vector
  void createRhsDiscretizationMatrix();
};

template<typename IntegratorType>
class FiniteElementMethodBaseTimeStepping<Mesh::RegularFixed<2ul>, BasisFunction::Lagrange<1>, IntegratorType> : 
  public FiniteElementMethodStiffnessMatrix<Mesh::RegularFixed<2ul>, BasisFunction::Lagrange<1>, IntegratorType>
{
public:
  using FiniteElementMethodStiffnessMatrix<Mesh::RegularFixed<2ul>, BasisFunction::Lagrange<1>, IntegratorType>::FiniteElementMethodStiffnessMatrix;

protected:
  //! create the discretization matrix which is the mapping between strong formulated and weak formulated rhs vector
  void createRhsDiscretizationMatrix();
};

template<typename IntegratorType>
class FiniteElementMethodBaseTimeStepping<Mesh::RegularFixed<3ul>, BasisFunction::Lagrange<1>, IntegratorType> : 
  public FiniteElementMethodStiffnessMatrix<Mesh::RegularFixed<3ul>, BasisFunction::Lagrange<1>, IntegratorType>
{
public:
  using FiniteElementMethodStiffnessMatrix<Mesh::RegularFixed<3ul>, BasisFunction::Lagrange<1>, IntegratorType>::FiniteElementMethodStiffnessMatrix;

protected:
  //! create the discretization matrix which is the mapping between strong formulated and weak formulated rhs vector
  void createRhsDiscretizationMatrix();
};

template<typename BasisFunctionType, typename IntegratorType>
class FiniteElementMethodBaseTimeStepping<Mesh::Deformable<1ul>, BasisFunctionType, IntegratorType> : 
  public FiniteElementMethodStiffnessMatrix<Mesh::Deformable<1ul>, BasisFunctionType, IntegratorType>
{
public:
  using FiniteElementMethodStiffnessMatrix<Mesh::Deformable<1ul>, BasisFunctionType, IntegratorType>::FiniteElementMethodStiffnessMatrix;

protected:
  //! create the discretization matrix which is the mapping between strong formulated and weak formulated rhs vector
  void createRhsDiscretizationMatrix();
};


template<typename BasisFunctionType, typename IntegratorType>
class FiniteElementMethodBaseTimeStepping<Mesh::Deformable<2ul>, BasisFunctionType, IntegratorType> : 
  public FiniteElementMethodStiffnessMatrix<Mesh::Deformable<2ul>, BasisFunctionType, IntegratorType>
{
public:
  using FiniteElementMethodStiffnessMatrix<Mesh::Deformable<2ul>, BasisFunctionType, IntegratorType>::FiniteElementMethodStiffnessMatrix;

protected:
  //! create the discretization matrix which is the mapping between strong formulated and weak formulated rhs vector
  void createRhsDiscretizationMatrix();
};

template<typename BasisFunctionType, typename IntegratorType>
class FiniteElementMethodBaseTimeStepping<Mesh::Deformable<3ul>, BasisFunctionType, IntegratorType> : 
  public FiniteElementMethodStiffnessMatrix<Mesh::Deformable<3ul>, BasisFunctionType, IntegratorType>
{
public:
  using FiniteElementMethodStiffnessMatrix<Mesh::Deformable<3ul>, BasisFunctionType, IntegratorType>::FiniteElementMethodStiffnessMatrix;

protected:
  //! create the discretization matrix which is the mapping between strong formulated and weak formulated rhs vector
  void createRhsDiscretizationMatrix();
};

// base class implementing timestepping as for diffusion equation
template<typename MeshType, typename BasisFunctionType, typename IntegratorType>
class FiniteElementMethodTimeStepping :
  public FiniteElementMethodBaseTimeStepping<MeshType, BasisFunctionType, IntegratorType>,
  public DiscretizableInTime
{
public:
  FiniteElementMethodTimeStepping(const DihuContext &context);
 
  //! proceed time stepping by computing output = stiffnessMatrix*input, output back in strong form
  void evaluateTimesteppingRightHandSide(Vec &input, Vec &output, int timeStepNo, double currentTime);
  
  //! initialize for use with timestepping
  void initialize() override;
  
  //! return true because the object has a specified mesh type
  bool knowsMeshType();
  
  friend class StiffnessMatrixTester;    ///< a class used for testing 
protected:
 
  //! do nothing, needed for initialize of base class that is overridden anyway
  void setRightHandSide(){};
  
  //! Extract from the rhs in weak formulation the rhs vector in strong formulation
  void recoverRightHandSide(Vec &result);
  
  double relativeTolerance_;      ///< relative tolerance for solver 
};



};  // namespace