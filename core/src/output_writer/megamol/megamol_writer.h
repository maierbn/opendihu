#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <vector>

#include "control/types.h"
#include "output_writer/generic.h"

namespace OutputWriter
{

template<typename FunctionSpaceType, typename OutputFieldVariablesType>
class MegaMolWriter
{
public:
  //! output the given field variable data for one mesh
  static void outputData(OutputFieldVariablesType fieldVariables, std::string meshName, std::shared_ptr<FunctionSpaceType> functionSpace,
                         PythonConfig specificSettings);
};

} // namespace

#include "output_writer/megamol/megamol_writer.tpp"
