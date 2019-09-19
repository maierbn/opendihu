#pragma once

#include <Python.h>  // has to be the first included header
namespace TimeSteppingScheme
{

/** Data type of output connector for QuasiStaticLinearElasticitySolver and others
 *  This is a separate class such that the mapping to the other solver object can be specialized.
 */
template<typename FieldVariableType>
struct ElasticitySolverOutputConnectorDataType
{
  std::shared_ptr<FieldVariableType> activation;   // field variable on 3D function space
};

//! output method for the ElasticitySolverOutputConnectorDataType type
template<typename FieldVariableType>
std::ostream &operator<<(std::ostream &stream, const ElasticitySolverOutputConnectorDataType<FieldVariableType> &rhs)
{
  stream << "<activation: " << *rhs.activation;
  return stream;
}


}  // namespace
