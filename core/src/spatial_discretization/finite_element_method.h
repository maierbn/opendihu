#pragma once

#include <Python.h>
#include <memory>

#include "spatial_discretization/spatial_discretization.h"
#include "time_stepping_scheme/discretizable_in_time.h"
#include "control/runnable.h"
#include "control/dihu_context.h"
#include "data_management/finite_elements.h"
#include "equation/laplace.h"
#include "equation/poisson.h"
#include "equation/type_traits.h"
#include "mesh/mesh.h"

#include "gtest/gtest_prod.h"

namespace SpatialDiscretization
{
 
template<typename MeshType, typename BasisFunctionType>
class FiniteElementMethodBase : public SpatialDiscretization, public Runnable
{
public:
  FiniteElementMethodBase(DihuContext &context);
  
  // perform computation
  void run();
  
  //! initialize for use with timestepping
  void initialize();
    
  //! get the stored mesh
  std::shared_ptr<Mesh::Mesh> mesh();
  
  friend class FunctionTester;    ///< a class used for testing 
  
protected:
 
  FRIEND_TEST(LaplaceTest, MatrixIsCorrect);
 
  virtual void setRightHandSide() = 0;
  virtual void applyBoundaryConditions();
  virtual void setStiffnessMatrix();
  virtual void solve();
  
  DihuContext &context_;    ///< the context object containing everything to be stored
  Data::FiniteElements data_;     ///< data object that holds all PETSc vectors and matrices
  PyObject *specificSettings_;    ///< python object containing the value of the python config dict with corresponding key
};
 
// inherited class that has additional Term template parameter
template<typename MeshType, typename BasisFunctionType, typename Term, typename DummyForTypeTraits = Term>
class FiniteElementMethod : public FiniteElementMethodBase<MeshType, BasisFunctionType>
{
};

// partial specialisation for Equation::Static::Laplace
template<typename MeshType, typename BasisFunctionType>
class FiniteElementMethod<MeshType, BasisFunctionType, Equation::Static::Laplace> :
  public FiniteElementMethodBase<MeshType, BasisFunctionType>
{
public:
  FiniteElementMethod(DihuContext &context);
 
private:
  void setRightHandSide();
};


// partial specialisation for Equation::hasLaplaceOperatorWithRhs, common base class for all meshtypes and basisfunction types
template<typename MeshType, typename BasisFunctionType>
class FiniteElementMethodBaseRhs :
  public FiniteElementMethodBase<MeshType, BasisFunctionType>, public DiscretizableInTime
{
public:
  FiniteElementMethodBaseRhs(DihuContext &context);
 
  //! Extract from the rhs in weak formulation the rhs vector in strong formulation
  void recoverRightHandSide(Vec &result);
  
  void evaluateTimesteppingRightHandSide(Vec &input, Vec &output);
  
  //! initialize for use with timestepping
  void initialize();
  
  //! get the stored mesh
  std::shared_ptr<Mesh::Mesh> mesh();
  
protected:
 
  //! Transform values in rhs vector into FEM discretized values by multiplying them with the integrate basis functions
  void multiplyRhsFactor();
  
  //! create the discretization matrix which is the mapping between strong formulated and weak formulated rhs vector
  void createDiscretizationMatrix();
 
  void setRightHandSide();
  
  bool timeSteppingInitialized_;    ///< if for use with timestepping the stiffness matrix and rhs vector are initialized
};

// common class for not specialized MeshType, BasisFunctionType
template<typename MeshType, typename BasisFunctionType, typename Term>
class FiniteElementMethod<MeshType, BasisFunctionType, Term, Equation::hasLaplaceOperatorWithRhs<Term>> :
  public FiniteElementMethodBaseRhs<MeshType, BasisFunctionType>
{
public:
  //! constructor
  FiniteElementMethod(DihuContext &context);
  
};

}  // namespace

#include "spatial_discretization/finite_element_method.tpp"
#include "spatial_discretization/finite_element_method_laplace.tpp"
#include "spatial_discretization/finite_element_method_poisson.tpp"