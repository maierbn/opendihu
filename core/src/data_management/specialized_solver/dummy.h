#pragma once

#include <Python.h>  // has to be the first included header

namespace Data
{

/**  The datastructures for the Dummy
  */
template<typename FunctionSpaceType>
class Dummy : public Data<FunctionSpaceType>
{
public:

  //! define the type of output connection variables, i.e. the values that will be transferred if the solver is part of a splitting or coupling scheme
  typedef OutputConnectorData<FunctionSpaceType,1,1> OutputConnectorDataType;

  //! constructor
  Dummy(DihuContext context);

  //! initialize and create all field variables
  void initialize();

  //! return the object that will be used to transfer values between solvers, in this case this includes only Vm
  std::shared_ptr<OutputConnectorDataType> getOutputConnectorData();

private:

  //! create all field variables with their respective sizes, this will be called automatically within initialize by the base class
  void createPetscObjects() override;

  std::shared_ptr<OutputConnectorDataType> outputConnectorData_;      //< the object that stores all components of field variables that will be transferred to other solvers
};

} // namespace Data

#include "data_management/specialized_solver/dummy.tpp"
