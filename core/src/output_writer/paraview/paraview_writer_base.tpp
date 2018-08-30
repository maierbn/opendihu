#include "output_writer/paraview/paraview_writer.h"

#include <fstream>
#include <iostream>
#include <thread>
#include <chrono>

#include "easylogging++.h"
#include "base64.h"

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "mesh/structured_regular_fixed.h"
#include "mesh/structured_deformable.h"
#include "mesh/unstructured_deformable.h"
#include "mesh/mesh.h"
#include "function_space/function_space.h"

namespace OutputWriter
{

template<typename DataType>
void ParaviewWriter::write(DataType& data, int timeStepNo, double currentTime)
{
  // check if output should be written in this timestep and prepare filename
  if (!Generic::prepareWrite(data, timeStepNo, currentTime))
  {
    return;
  }

  writeSolutionDim<DataType::FunctionSpace::dim()>(data);
}

template <int dimension, typename DataType>
void ParaviewWriter::writeSolutionDim(DataType &data)
{
  LOG(TRACE) << "writeMesh<" <<dimension<< ">()";

  if (std::dynamic_pointer_cast<Mesh::StructuredRegularFixedOfDimension<dimension>>(data.functionSpace()) != NULL)
  {
    writeRectilinearGrid<Mesh::StructuredRegularFixedOfDimension<dimension>>(data);
  }
  else if (std::dynamic_pointer_cast<
             FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<dimension>,BasisFunction::LagrangeOfOrder<1>>
           >(data.functionSpace()) != NULL)
  {
    // structured grid only for elements that contain only nodes at the corners (i.e. linear lagrange elements)
    writeStructuredGrid<dimension>(data);
  }
  else if (std::dynamic_pointer_cast<Mesh::UnstructuredDeformableOfDimension<dimension>>(data.functionSpace()) != NULL)
  {
    writeUnstructuredGrid<dimension>(data);
  }
}

};
