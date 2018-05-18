#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <vector>

namespace OutputWriter
{

/** This is a base class for python output writer and contains functionality to create a python dictionary from some field variables
 * that will then be output.
 */
template<typename OutputFieldVariablesType>
class PythonBase
{
public:
  //! create a python dict that contains data and meta data of field variables
  //! @param onlyNodalValues: if only values at nodes should be contained, this discards the derivative values for Hermite
  static PyObject *buildPyFieldVariablesObject(OutputFieldVariablesType fieldVariables, bool onlyNodalValues);

  //! create a python object for the given field variable
  template<typename FieldVariableType>
  static PyObject *buildPyFieldVariableObject(FieldVariableType fieldVariable, bool onlyNodalValues);
};

namespace PythonLoopOverTuple
{

 /** Static recursive loop from 0 to number of entries in the tuple
 *  Stopping criterion
 */
template<typename OutputFieldVariablesType, int i=0>
inline typename std::enable_if<i == std::tuple_size<OutputFieldVariablesType>::value, void>::type
loopOverTuple(const OutputFieldVariablesType &fieldVariables, PyObject *pyData, bool onlyNodalValues)
{}

/** Static recursive loop from 0 to number of entries in the tuple
 *  Loop body
 */
template<typename OutputFieldVariablesType, int i=0>
inline typename std::enable_if<i < std::tuple_size<OutputFieldVariablesType>::value, void>::type
loopOverTuple(const OutputFieldVariablesType &fieldVariables, PyObject *pyData, bool onlyNodalValues);

};   // namespace

};  // namespace

#include "output_writer/python/python_base.tpp"
