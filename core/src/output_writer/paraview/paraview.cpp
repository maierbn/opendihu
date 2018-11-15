#include "output_writer/paraview/paraview.h"

#include <fstream>
#include <iostream>
#include <thread>
#include <chrono>

#include "easylogging++.h"
#include "base64.h"

#include <utility/python_utility.h>
#include <utility/petsc_utility.h>
#include <mesh/structured_regular_fixed.h>
#include <mesh/structured_deformable.h>
#include <mesh/unstructured_deformable.h>
#include <mesh/mesh.h>

namespace OutputWriter
{

Paraview::Paraview(DihuContext context, PythonConfig settings) :
  Generic(context, settings)
{
  binaryOutput_ = settings.getOptionBool("binary", true);
  fixedFormat_ = settings.getOptionBool("fixedFormat", true);
  combineFiles_ = settings.getOptionBool("combineFiles", false);
}

std::string Paraview::encodeBase64Vec(const Vec &vector, bool withEncodedSizePrefix)
{
  int vectorSize = 0;
  VecGetSize(vector, &vectorSize);

  std::vector<int> indices(vectorSize);
  std::iota (indices.begin(), indices.end(), 0);    // fill with increasing numbers: 0,1,2,...
  std::vector<double> values(vectorSize);
  VecGetValues(vector, vectorSize, indices.data(), values.data());

  return encodeBase64Float(values.begin(), values.end(), withEncodedSizePrefix);
}

std::string Paraview::convertToAscii(const Vec &vector, bool fixedFormat)
{
  std::vector<double> vectorValues;
  PetscUtility::getVectorEntries(vector, vectorValues);

  return convertToAscii(vectorValues, fixedFormat);
}

std::string Paraview::convertToAscii(const std::vector<double> &vector, bool fixedFormat)
{
  std::stringstream result;
  for (auto value : vector)
  {
    if (fixedFormat)
    {
      result << std::setw(16) << std::scientific << value << " ";
    }
    else
    {
      result << value << " ";
    }
  }
  return result.str();
}

std::string Paraview::convertToAscii(const std::vector<int> &vector, bool fixedFormat)
{
  std::stringstream result;
  for (auto value : vector)
  {
    if (fixedFormat)
    {
      result << std::setw(16) << std::scientific << (float)(value) << " ";
    }
    else
    {
      result << value << " ";
    }
  }
  return result.str();
}
};
