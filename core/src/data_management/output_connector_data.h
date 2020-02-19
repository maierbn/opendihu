#pragma once

#include <Python.h>  // has to be the first included header
#include "field_variable/field_variable.h"

namespace Data
{

/** The data type for the outputConnector of the TimeStepping class.
  *  This is the data that will be transferred to connected classes in the nested structure.
  */
template<typename FunctionSpaceType, int nComponents>
struct ComponentOfFieldVariable
{
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> values;    //< a field variable containing the payload data that is to be exchangend to another solver
  int componentNo;                              //< the component of values that is relevant, only this component out of the potentially multi-component field variable in values will be transferred.

  //! constructor
  ComponentOfFieldVariable(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> fieldVariable, int componentNo);

  //! setValue convenience method, sets the value at the given dofNoLocal
  void setValue(dof_no_t dofNoLocal, double value, InsertMode petscInsertMode=INSERT_VALUES);

  //! set all values
  void setValuesWithoutGhosts(const std::vector<double> &values, InsertMode petscInsertMode=INSERT_VALUES);

  //! get the value of the local dof no
  double getValue(dof_no_t dofNoLocal);


};

// operator used for output
template<typename FunctionSpaceType, int nComponents>
std::ostream &operator<<(std::ostream &stream, const ComponentOfFieldVariable<FunctionSpaceType,nComponents> &rhs);

/** This is a general data type that can contain multiple field variables with designated components.
 *  They are used for the outputConnector of different solvers.
 *  This particular struct is, for example, used for the CellMLAdapter, where the two variables store states and intermediates to be passed to the next solver.
 */
template<typename FunctionSpaceType, int nComponents1, int nComponents2=1>
struct OutputConnectorData
{
  std::vector<ComponentOfFieldVariable<FunctionSpaceType,nComponents1>> variable1;    //< vector of indications of components of field variables, the field variables have all the same number of components
  std::vector<ComponentOfFieldVariable<FunctionSpaceType,nComponents2>> variable2;    //< second vector with different number of components for the field variables

  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> geometryField = nullptr;   //< geometry field to be transferred if set

  //! add a component of a field variable to the vector
  void addFieldVariable(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents1>> fieldVariable, int componentNo=0);

  //! add a component of a field variable to the vector
  void addFieldVariable2(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents2>> fieldVariable, int componentNo=0);

  //! assign a geometry field, this means that geometry data will be transferred
  void addGeometryField(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> geometryField);
};

// operator used for output
template<typename FunctionSpaceType, int nComponents1, int nComponents2>
std::ostream &operator<<(std::ostream &stream, const OutputConnectorData<FunctionSpaceType,nComponents1,nComponents2> &rhs);

// operator used for output
template<typename FunctionSpaceType, int nComponents1, int nComponents2>
std::ostream &operator<<(std::ostream &stream, const std::shared_ptr<OutputConnectorData<FunctionSpaceType,nComponents1,nComponents2>> &rhs);

} // namespace

#include "data_management/output_connector_data.tpp"
