#pragma once

#include <Python.h>  // has to be the first included header
#include "field_variable/field_variable.h"

namespace Data
{

/** The data type for the slotConnector of the TimeStepping class.
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
 *  They are used for the slotConnector of different solvers.
 *  This particular struct is, for example, used for the CellMLAdapter, where the two variables store states and algebraics to be passed to the next solver.
 */
template<typename FunctionSpaceType, int nComponents1, int nComponents2=1>
struct SlotConnectorData
{
  std::vector<ComponentOfFieldVariable<FunctionSpaceType,nComponents1>> variable1;    //< vector of indications of components of field variables, the field variables have all the same number of components
  std::vector<ComponentOfFieldVariable<FunctionSpaceType,nComponents2>> variable2;    //< second vector with different number of components for the field variables

  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> geometryField = nullptr;   //< geometry field to be transferred if set

  std::vector<std::string> slotNames;                                                 //< names for all slots of variable1 and variable2, used for connecting slots by their name, not number

  //! define types of the two possible field variables
  typedef FieldVariable::FieldVariable<FunctionSpaceType,nComponents1> FieldVariable1Type;
  typedef FieldVariable::FieldVariable<FunctionSpaceType,nComponents2> FieldVariable2Type;

  //! add a component of a field variable to the vector
  void addFieldVariable(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents1>> fieldVariable, int componentNo=0);

  //! add a component of a field variable to the vector
  void addFieldVariable2(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents2>> fieldVariable, int componentNo=0);

  //! assign a geometry field, this means that geometry data will be transferred
  void addGeometryField(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> geometryField);

  //! get the number of slots that is contained in this SlotConnectorData, this is equal to variable1.size() + variable2.size()
  int nSlots();
};

// operator used for output
template<typename FunctionSpaceType, int nComponents1, int nComponents2>
std::ostream &operator<<(std::ostream &stream, const SlotConnectorData<FunctionSpaceType,nComponents1,nComponents2> &rhs);

// operator used for output
template<typename FunctionSpaceType, int nComponents1, int nComponents2>
std::ostream &operator<<(std::ostream &stream, const std::shared_ptr<SlotConnectorData<FunctionSpaceType,nComponents1,nComponents2>> &rhs);

} // namespace

#include "slot_connection/slot_connector_data.tpp"
