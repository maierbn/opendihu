#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <vector>

#include "control/types.h"
#include "output_writer/generic.h"

namespace OutputWriter
{

/** Base class of ParaviewWriter that writes vtk files of given field variables.
 *  OutputFieldVariablesType is a std::tuple<std::shared_ptr<>, std::shared_ptr<>, ...> of field variables.
 *  Only field variables which are defined on the specified mesh will be output.
 *  The FunctionSpaceType has to be the type of the field variables given in meshName.
 */
template<typename FunctionSpaceType, typename OutputFieldVariablesType>
class ParaviewWriter
{
public:
  //! write paraview file to given filename, only output fieldVariables that are on a mesh with the given meshName 
  static void outputFile(std::string filename, OutputFieldVariablesType fieldVariables, std::string meshName, std::shared_ptr<FunctionSpaceType> mesh, int nFieldVariablesOfMesh);
  
private:

  //! write out solution templated by dimension
  template <int dimension, typename DataType>
  static void writeSolutionDim(DataType &data);
};

} // namespace

#include "output_writer/paraview/paraview_writer.tpp"
