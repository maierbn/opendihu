#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <vector>

#include "control/types.h"
#include "output_writer/generic.h"
#include "output_writer/megamol/megamol.h"

namespace OutputWriter
{

template<typename FunctionSpaceType, typename OutputFieldVariablesType>
class MegaMolWriter
{
public:

#ifdef HAVE_ADIOS
  //! output the given field variable data for one mesh
  static void outputData(OutputFieldVariablesType fieldVariables, std::string meshName, std::shared_ptr<FunctionSpaceType> functionSpace,
                         PythonConfig specificSettings, std::shared_ptr<adios2::Engine> adiosWriter, std::shared_ptr<adios2::IO> adiosIo, BoundingBox &boundingBox);
private:

  static std::map<std::string, adios2::Variable<double>> geometryTable_;    // the geometry data for each mesh
  static std::map<std::string, adios2::Variable<double>> vmTable_;          // the vm datat for each mesh

  static std::map<std::string, std::shared_ptr<std::vector<double>>> geometryTableData_;   ///< buffer for the values to be written by ADIOS
  static std::map<std::string, std::shared_ptr<std::vector<double>>> vmTableData_;   ///< buffer for the values to be written by ADIOS
#endif
};

} // namespace

#include "output_writer/megamol/megamol_writer.tpp"
