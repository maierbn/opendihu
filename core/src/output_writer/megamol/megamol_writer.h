#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <vector>

#include "control/types.h"
#include "output_writer/generic.h"
#include "output_writer/megamol/megamol.h"

namespace OutputWriter
{

template<typename FunctionSpaceType, typename FieldVariablesForOutputWriterType>
class MegaMolWriter
{
public:

#ifdef HAVE_ADIOS
  //! output the given field variable data for one mesh
  static void outputData(FieldVariablesForOutputWriterType fieldVariables, std::string meshName, std::shared_ptr<FunctionSpaceType> functionSpace,
                         PythonConfig specificSettings, MegaMolWriterContext &megaMolWriterContext);
#endif
};

} // namespace

#include "output_writer/megamol/megamol_writer.tpp"
