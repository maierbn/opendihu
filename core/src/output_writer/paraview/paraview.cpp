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

Paraview::Paraview(PyObject *settings) : Generic(settings)
{
}

std::string Paraview::encodeBase64(const Vec &vector)
{
  int vectorSize = 0;
  VecGetSize(vector, &vectorSize);

  std::vector<int> indices(vectorSize);
  std::iota (indices.begin(), indices.end(), 0);    // fill with increasing numbers: 0,1,2,...
  std::vector<double> values(vectorSize);
  VecGetValues(vector, vectorSize, indices.data(), values.data());

  return encodeBase64(values);
}

std::string Paraview::convertToAscii(const Vec &vector, bool fixedFormat)
{
  std::vector<double> vectorValues;
  PetscUtility::getVectorEntries(vector, vectorValues);

  return convertToAscii(vectorValues, fixedFormat);
}

std::string Paraview::encodeBase64(const std::vector<double> &vector)
{
  // encode as Paraview Float32
  assert(sizeof(float) == 4);
  
  int rawLength = vector.size()*sizeof(float);
  int encodedLength = Base64::EncodedLength(4+rawLength);

  char raw[4+rawLength];
  for (unsigned int i=0; i<vector.size(); i++)
  {
    union {
      float d;
      char c[4];
    };
    d = vector[i];
    memcpy(raw+4+i*sizeof(float), c, 4);
  }

  // prepend number of bytes as uint32
  union {
    uint32_t i;
    char c[4];
  };
  i = rawLength;

  memcpy(raw, c, 4);

  char encoded[encodedLength+1];
  //Base64::Encode(reinterpret_cast<char *>(vector.data()), rawLength, encoded, encodedLength);
  bool success = Base64::Encode(raw, rawLength+4, encoded, encodedLength);
  if (!success)
    LOG(WARNING) << "encoding failed";

  encoded[encodedLength] = '\0';

  return std::string(encoded);
}

std::string Paraview::encodeBase64(const std::vector<int> &vector)
{
  // encode as Paraview Int32
  assert(sizeof(int) == 4);

  int rawLength = vector.size()*sizeof(int);
  int encodedLength = Base64::EncodedLength(4+rawLength);

  char raw[4+rawLength];
  for (unsigned int i=0; i<vector.size(); i++)
  {
    union {
      int integer;
      char c[4];
    };
    integer = vector[i];
    memcpy(raw+4+i*sizeof(int), c, 4);
  }

  // prepend number of bytes as uint32
  union {
    uint32_t i;
    char c[4];
  };
  i = rawLength;

  memcpy(raw, c, 4);

  char encoded[encodedLength+1];
  //Base64::Encode(reinterpret_cast<char *>(vector.data()), rawLength, encoded, encodedLength);
  bool success = Base64::Encode(raw, rawLength+4, encoded, encodedLength);
  if (!success)
    LOG(WARNING) << "encoding failed";

  encoded[encodedLength] = '\0';

  return std::string(encoded);
}

std::string Paraview::convertToAscii(const std::vector<double> &vector, bool fixedFormat)
{
  std::stringstream result;
  for(auto value : vector)
  {
    if(fixedFormat)
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
  for(auto value : vector)
  {
    if(fixedFormat)
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
