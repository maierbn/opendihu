#pragma once

#include <Python.h>  // has to be the first included header
#include <petscsys.h>
#include <petscksp.h>
#include <iostream>
#include <vector>

#include "data_management/solution_vector_mapping.h"

class DiscretizableInTime
{
public:
  //! constructor
  DiscretizableInTime();

  // Classes that derive from DiscretizableInTime must define a constexpr nComponents that specifies the number of components in the solution field variable
  //typedef .. nComponents;

  //! initialize timestepping
  virtual void initialize() = 0;
  
  //! initialize timestepping
  virtual void initializeForImplicitTimeStepping() = 0;

  //! timestepping rhs function f of equation u_t = f(u,t), compute output = f(u,t) where u=input
  virtual void evaluateTimesteppingRightHandSideExplicit(Vec &input, Vec &output, int timeStepNo, double currentTime) = 0;
  
  //! timestepping rhs function f of equation u_t = f(u,t), computed in an implicit time step with equation Au^{t+1}=f^{t} where A=(I-dt*M^{-1}K). This is called when implicit Euler is used.
  //virtual void evaluateTimesteppingRightHandSideImplicit(Vec &input, Vec &output, int timeStepNo, double currentTime) = 0;

  //! get the number of degrees of freedom per node which is 1 by default
  virtual int nComponentsNode();

  //! set initial values and return true or don't do anything and return false
  // this could use std::any (c++17)
  //template<typename FieldVariableType>
  //virtual bool setInitialValues(std::shared_ptr<FieldVariableType> initialValues);

  //! get the names of components to be used for the solution variable
  virtual void getComponentNames(std::vector<std::string> &componentNames);
  
  //! set the subset of ranks that will compute the work
  virtual void setRankSubset(Partition::RankSubset rankSubset) = 0;
  
  //! set if the class should handle dirichlet boundary conditions. A time stepping scheme sets this to false, because for dynamic problems the time stepping scheme handles the boundary conditions, not e.g. the FiniteElementMethod.
  //! By default it is set to true, which is needed for static problems, like Laplace.
  virtual void setBoundaryConditionHandlingEnabled(bool boundaryConditionHandlingEnabled) = 0;

  //! return whether the object has a specified mesh type and is not independent of the mesh type
  virtual bool knowsMeshType() = 0;

  //! return the mesh
  //virtual std::shared_ptr<FunctionSpaceType> functionSpace() = 0;
  //old: virtual std::shared_ptr<Mesh::Mesh> mesh() = 0;

protected:
};
