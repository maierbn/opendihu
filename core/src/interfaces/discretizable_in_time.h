#pragma once

#include <Python.h>  // has to be the first included header
#include <petscsys.h>
#include <petscksp.h>
#include <iostream>
#include <vector>

#include "partition/mesh_partition/01_mesh_partition.h"

template<typename FunctionSpaceType, int nComponentsSolutionVariable>
class DiscretizableInTime
{
public:
  // Classes that derive from DiscretizableInTime must define a constexpr nComponents that specifies the number of components in the solution field variable
  //! get the number of components
  static constexpr int nComponents();

  //! initialize timestepping
  virtual void initialize() = 0;
  
  //! initialize timestepping
  virtual void initializeForImplicitTimeStepping() = 0;

  //! timestepping rhs function f of equation u_t = f(u,t), compute output = f(u,t) where u=input
  virtual void evaluateTimesteppingRightHandSideExplicit(Vec &input, Vec &output, int timeStepNo, double currentTime) = 0;
  
  //! timestepping rhs function f of equation u_t = f(u,t), computed in an implicit time step with equation Au^{t+1}=f^{t} where A=(I-dt*M^{-1}K). This is called when implicit Euler is used.
  //virtual void evaluateTimesteppingRightHandSideImplicit(Vec &input, Vec &output, int timeStepNo, double currentTime) = 0;

  //! Get the data that will be transferred in the operator splitting to the other term of the splitting.
  //! The transfer is done by the slot_connector_data_transfer class.
  //virtual std::shared_ptr<SlotConnectorDataType> getSlotConnectorData();

  //! this will be called right before getSlotConnectorData
  virtual void prepareForGetSlotConnectorData();

  //! set initial values and return true or don't do anything and return false
  virtual bool setInitialValues(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponentsSolutionVariable>> initialValues) = 0;

  //! get the names of components to be used for the solution variable
  virtual void getComponentNames(std::vector<std::string> &componentNames);
  
  //! set the subset of ranks that will compute the work
  virtual void setRankSubset(Partition::RankSubset rankSubset) = 0;
  
  //! set if the class should handle dirichlet boundary conditions. A time stepping scheme sets this to false, because for dynamic problems the time stepping scheme handles the boundary conditions, not e.g. the FiniteElementMethod.
  //! By default it is set to true, which is needed for static problems, like Laplace.
  virtual void setBoundaryConditionHandlingEnabled(bool boundaryConditionHandlingEnabled) = 0;

   //! set the solution field variable in the data object, that actual data is stored in the timestepping scheme object
  virtual void setSolutionVariable(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponentsSolutionVariable>> solution) = 0;

  //! pass on the slot connector data object from the timestepping scheme object to be modified,
  //! this is needed for other DiscretizableInTime objects
  virtual void setSlotConnectorData(std::shared_ptr<Data::SlotConnectorData<FunctionSpaceType,nComponentsSolutionVariable>> slotConnectorDataTimeStepping) = 0;


protected:
};

#include "interfaces/discretizable_in_time.tpp"
