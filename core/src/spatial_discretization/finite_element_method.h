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
#include "output_writer/manager.h"

namespace SpatialDiscretization
{
 
template<typename MeshType, typename BasisFunctionType>
class FiniteElementMethodBase : public SpatialDiscretization, public Runnable
{
public:
  FiniteElementMethodBase(const DihuContext &context);
  
  // perform computation
  void run();
  
  //! initialize for use as laplace or poisson equation, not for timestepping
  void initialize();
    
  //! get the stored mesh
  std::shared_ptr<Mesh::Mesh> mesh();
  
  friend class StiffnessMatrixTester;    ///< a class used for testing 
protected:
 
  virtual void setRightHandSide() = 0;
  virtual void applyBoundaryConditions();
  virtual void setStiffnessMatrix();
  virtual void solve();
  
  const DihuContext &context_;    ///< the context object containing everything to be stored
  Data::FiniteElements data_;     ///< data object that holds all PETSc vectors and matrices
  PyObject *specificSettings_;    ///< python object containing the value of the python config dict with corresponding key
  OutputWriter::Manager outputWriterManager_; ///< manager object holding all output writer
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
  FiniteElementMethod(const DihuContext &context);
 
private:
  void setRightHandSide();
};


// base class implementing right hand side, that can be set by user for poisson equation
template<typename MeshType, typename BasisFunctionType>
class FiniteElementMethodBaseRhs :
  public FiniteElementMethodBase<MeshType, BasisFunctionType>
{
public:
  FiniteElementMethodBaseRhs(const DihuContext &context);
 
  friend class StiffnessMatrixTester;    ///< a class used for testing 
protected:
 
  //! Transform values in rhs vector into FEM discretized values by multiplying them with the integrate basis functions
  void transferRhsToWeakForm();
  
  //! read in rhs values from config and transfer to weak form
  void setRightHandSide();
};

// base class implementing timestepping as for diffusion equation
template<typename MeshType, typename BasisFunctionType>
class FiniteElementMethodBaseTimeStepping :
  public FiniteElementMethodBase<MeshType, BasisFunctionType>, 
  public DiscretizableInTime
{
public:
  FiniteElementMethodBaseTimeStepping(const DihuContext &context);
 
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
  
  //! create the discretization matrix which is the mapping between strong formulated and weak formulated rhs vector
  void createRhsDiscretizationMatrix();
  
  double relativeTolerance_;      ///< relative tolerance for solver 
};



// common class for not specialized MeshType, BasisFunctionType, for poisson equation
template<typename MeshType, typename BasisFunctionType, typename Term>
class FiniteElementMethod<MeshType, BasisFunctionType, Term, Equation::hasLaplaceOperatorWithRhs<Term>> :
  public FiniteElementMethodBaseRhs<MeshType, BasisFunctionType>
{
public:
  //! constructor
  FiniteElementMethod(const DihuContext &context);
  
};

// common class for not specialized MeshType, BasisFunctionType, for time stepping
template<typename MeshType, typename BasisFunctionType, typename Term>
class FiniteElementMethod<MeshType, BasisFunctionType, Term, Equation::hasLaplaceOperatorWithTimeStepping<Term>> :
  public FiniteElementMethodBaseTimeStepping<MeshType, BasisFunctionType>
{
public:
  //! constructor
  FiniteElementMethod(const DihuContext &context);
  
};


}  // namespace

#include "spatial_discretization/finite_element_method.tpp"
#include "spatial_discretization/finite_element_method_laplace.tpp"
#include "spatial_discretization/finite_element_method_poisson.tpp"
#include "spatial_discretization/finite_element_method_timestepping.tpp"