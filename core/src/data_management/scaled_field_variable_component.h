#pragma once

#include <Python.h>  // has to be the first included header
#include "field_variable/field_variable.h"

namespace Data
{


/** The data type for the outputConnector of the TimeStepping class.
  *  This is the data that will be transferred to connected classes in the nested structure.
  */
template<typename FunctionSpaceType, int nComponents>
struct ScaledFieldVariableComponent
{
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> values;    //< a field variable containing the payload data that is to be exchangend to another solver
  int componentNo;                              //< the component of values that is relevant, only this component out of the potentially multi-component field variable in values will be transferred.
  double scalingFactor;    //< a scaling factor, the values will be multiplied by the factor before the transfer. Disabled if set to 1.0.
};

// operator used for output
template<typename FunctionSpaceType, int nComponents>
std::ostream &operator<<(std::ostream &stream, const ScaledFieldVariableComponent<FunctionSpaceType,nComponents> &rhs)
{
  stream << "<" << *(rhs.values) << " (componentNo " << rhs.componentNo << ", scalingFactor " << rhs.scalingFactor << ")>";
  return stream;
}

} // namespace
