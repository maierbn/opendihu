#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "control/types.h"
#include "mesh/mesh.h"
#include "data_management/data.h"
#include "data_management/time_stepping/time_stepping.h"
#include "field_variable/field_variable.h"

namespace Data
{

/**  The datastructures for the PrescribedValues dummy
  */
template<typename FunctionSpaceType, int nComponents1=1, int nComponents2=1>
class PrescribedValues : public Data<FunctionSpaceType>
{
public:

  //! define the type of the first field variable for convenience
  typedef FieldVariable::FieldVariable<FunctionSpaceType,nComponents1> FieldVariable1Type;

  //! define the type of the second field variable for convenience
  typedef FieldVariable::FieldVariable<FunctionSpaceType,nComponents2> FieldVariable2Type;

  //! define the type of output connection variables, i.e. the values that will be transferred if the solver is part of a splitting or coupling scheme
  typedef OutputConnectorData<FunctionSpaceType,nComponents1,nComponents2> OutputConnectorDataType;

  //! constructor
  PrescribedValues(DihuContext context);

  //! return a reference to the fieldVariable with given index with nComponents1 components
  std::shared_ptr<FieldVariable1Type> fieldVariable1(int index);

  //! return a reference to the fieldVariable with given index with nComponents2 components
  std::shared_ptr<FieldVariable2Type> fieldVariable2(int index);

  //! initialize and create all field variables
  void initialize(std::vector<std::string> fieldVariable1Names, std::vector<std::string> fieldVariable2Names);

  //! print all stored data to stdout
  void print();

  //! return the object that will be used to transfer values between solvers, in this case this includes only Vm
  std::shared_ptr<OutputConnectorDataType> getOutputConnectorData();

  //! field variables that will be output by outputWriters
  typedef std::tuple<
    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>>,     // geometry, this always has to be the first field variable, such that the output writer knows the geometry of the mesh
    std::vector<std::shared_ptr<FieldVariable1Type>>,     // variable1,
    std::vector<std::shared_ptr<FieldVariable2Type>>      // variable2
    // ... add all field variables that you want to have in the output file
  > FieldVariablesForOutputWriter;

  //! get pointers to all field variables that can be written by output writers
  FieldVariablesForOutputWriter getFieldVariablesForOutputWriter();

private:

  //! create all field variables with their respective sizes, this will be called automatically within initialize by the base class
  void createPetscObjects() override;

  std::vector<std::shared_ptr<FieldVariable1Type>> fieldVariables1_;   //< the first list of field variables, with nComponents1 components
  std::vector<std::shared_ptr<FieldVariable2Type>> fieldVariables2_;   //< the second list of field variables, with nComponents2 components

  std::shared_ptr<OutputConnectorDataType> outputConnectorData_;      //< the object that stores all components of field variables that will be transferred to other solvers

  std::vector<std::string> fieldVariable1Names_;                      //< names of the field variables with `nComponents1` components
  std::vector<std::string> fieldVariable2Names_;                      //< names of the field variables with `nComponents1` components

};

} // namespace Data

#include "data_management/specialized_solver/prescribed_values.tpp"
