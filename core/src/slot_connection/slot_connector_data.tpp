#include "slot_connection/slot_connector_data.h"

namespace Data
{

template<typename FunctionSpaceType, int nComponents>
ComponentOfFieldVariable<FunctionSpaceType,nComponents>::
ComponentOfFieldVariable(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> fieldVariable, int _componentNo)
{
  values = fieldVariable;
  componentNo = _componentNo;
}

template<typename FunctionSpaceType, int nComponents>
void ComponentOfFieldVariable<FunctionSpaceType,nComponents>::
setValue(dof_no_t dofNoLocal, double value, InsertMode petscInsertMode)
{
  assert(values->nComponents() > componentNo);
  values->setValue(componentNo, dofNoLocal, value, petscInsertMode);
}

template<typename FunctionSpaceType, int nComponents>
void ComponentOfFieldVariable<FunctionSpaceType,nComponents>::
setValuesWithoutGhosts(const std::vector<double> &values, InsertMode petscInsertMode)
{
  assert(this->values->nComponents() > componentNo);
  this->values->setValuesWithoutGhosts(componentNo, values, petscInsertMode);
}

template<typename FunctionSpaceType, int nComponents>
double ComponentOfFieldVariable<FunctionSpaceType,nComponents>::
getValue(dof_no_t dofNoLocal)
{
  assert(values->nComponents() > componentNo);
  return values->getValue(componentNo, dofNoLocal);
}

// operator used for output
template<typename FunctionSpaceType, int nComponents>
std::ostream &operator<<(std::ostream &stream, const ComponentOfFieldVariable<FunctionSpaceType,nComponents> &rhs)
{
  stream << "(" << *(rhs.values) << " (componentNo " << rhs.componentNo << "))";
  return stream;
}

template<typename FunctionSpaceType, int nComponents1, int nComponents2>
void SlotConnectorData<FunctionSpaceType,nComponents1,nComponents2>::
addFieldVariable(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents1>> fieldVariable, int componentNo)
{
  variable1.push_back(ComponentOfFieldVariable<FunctionSpaceType,nComponents1>(fieldVariable, componentNo));
}

template<typename FunctionSpaceType, int nComponents1, int nComponents2>
void SlotConnectorData<FunctionSpaceType,nComponents1,nComponents2>::
addFieldVariable2(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents2>> fieldVariable, int componentNo)
{
  variable2.push_back(ComponentOfFieldVariable<FunctionSpaceType,nComponents2>(fieldVariable, componentNo));
}

template<typename FunctionSpaceType, int nComponents1, int nComponents2>
void SlotConnectorData<FunctionSpaceType,nComponents1,nComponents2>::
addGeometryField(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> _geometryField)
{
  geometryField = _geometryField;
}

template<typename FunctionSpaceType, int nComponents1, int nComponents2>
int SlotConnectorData<FunctionSpaceType,nComponents1,nComponents2>::
nSlots()
{
  return variable1.size() + variable2.size();
}

// operator used for output
template<typename FunctionSpaceType, int nComponents1, int nComponents2>
std::ostream &operator<<(std::ostream &stream, const SlotConnectorData<FunctionSpaceType,nComponents1,nComponents2> &rhs)
{
  if (rhs.variable1.empty())
    stream << "\t(variable1 empty\n";
  else
  {
    stream << "\t(variable1:\n";
    for (const ComponentOfFieldVariable<FunctionSpaceType,nComponents1> &entry : rhs.variable1)
    {
      stream << "\t\t[" << entry.values << ": " << *(entry.values) << " name \"" << entry.values->name() << "\", component " << entry.componentNo << "], " << std::endl;
    }
  }

  if (rhs.variable2.empty())
    stream << "\t, variable2 empty\n";
  else
  {
    stream << "\t, variable2:\n";
    for (const ComponentOfFieldVariable<FunctionSpaceType,nComponents2> &entry : rhs.variable2)
    {
      stream << "\t\t[" << entry.values << ": " << *(entry.values) << " name \"" << entry.values->name() << "\", component " << entry.componentNo << "], " << std::endl;
    }
  }
  stream << "\t)";
  return stream;
}

// operator used for output
template<typename FunctionSpaceType, int nComponents1, int nComponents2>
std::ostream &operator<<(std::ostream &stream, const std::shared_ptr<SlotConnectorData<FunctionSpaceType,nComponents1,nComponents2>> &rhs)
{
  //stream << rhs.get();
  if (rhs)
    stream << *rhs;
  else
  {
    // the program should not have any SlotConnectorData pointer uninitialized any time after initialize()
    LOG(FATAL) << "SlotConnectorData is null!";
  }
  return stream;
}

} // namespace
