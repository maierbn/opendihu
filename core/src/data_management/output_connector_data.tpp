#include "data_management/output_connector_data.h"

namespace Data
{

template<typename FunctionSpaceType, int nComponents>
ComponentOfFieldVariable<FunctionSpaceType,nComponents>::
ComponentOfFieldVariable(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> fieldVariable, int _componentNo)
{
  values = fieldVariable;
  componentNo = _componentNo;
}

// operator used for output
template<typename FunctionSpaceType, int nComponents>
std::ostream &operator<<(std::ostream &stream, const ComponentOfFieldVariable<FunctionSpaceType,nComponents> &rhs)
{
  stream << "<" << *(rhs.values) << " (componentNo " << rhs.componentNo << ", scalingFactor " << rhs.scalingFactor << ")>";
  return stream;
}

template<typename FunctionSpaceType, int nComponents1, int nComponents2>
void OutputConnectorData<FunctionSpaceType,nComponents1,nComponents2>::
addFieldVariable(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents1>> fieldVariable, int componentNo)
{
  variable1.push_back(ComponentOfFieldVariable<FunctionSpaceType,nComponents1>(fieldVariable, componentNo));
}

template<typename FunctionSpaceType, int nComponents1, int nComponents2>
void OutputConnectorData<FunctionSpaceType,nComponents1,nComponents2>::
addFieldVariable2(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents2>> fieldVariable, int componentNo)
{
  variable2.push_back(ComponentOfFieldVariable<FunctionSpaceType,nComponents2>(fieldVariable, componentNo));
}

template<typename FunctionSpaceType, int nComponents1, int nComponents2>
void OutputConnectorData<FunctionSpaceType,nComponents1,nComponents2>::
addGeometryField(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> _geometryField)
{
  geometryField = _geometryField;
}

// operator used for output
template<typename FunctionSpaceType, int nComponents1, int nComponents2>
std::ostream &operator<<(std::ostream &stream, const OutputConnectorData<FunctionSpaceType,nComponents1,nComponents2> &rhs)
{
  stream << "<";
  for (const ComponentOfFieldVariable<FunctionSpaceType,nComponents1> &entry : rhs.variable1)
  {
    stream << "[" << *(entry.values) << " component " << entry.componentNo << "], ";
  }
  stream << ";";
  for (const ComponentOfFieldVariable<FunctionSpaceType,nComponents2> &entry : rhs.variable2)
  {
    stream << "[" << *(entry.values) << " component " << entry.componentNo << "], ";
  }
  return stream;
}

} // namespace
