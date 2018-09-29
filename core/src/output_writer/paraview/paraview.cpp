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

Paraview::Paraview(DihuContext context, PyObject *settings) : Generic(context, settings)
{
  binaryOutput_ = PythonUtility::getOptionBool(settings, "binary", true);
  fixedFormat_ = PythonUtility::getOptionBool(settings, "fixedFormat", true);
  combineFiles_ = PythonUtility::getOptionBool(settings, "combineFiles", false);
}

void Paraview::setRankSubset(Partition::RankSubset &rankSubset)
{
  rankSubset_ = rankSubset;
}

std::string Paraview::encodeBase64(const Vec &vector, bool withEncodedSizePrefix)
{
  int vectorSize = 0;
  VecGetSize(vector, &vectorSize);

  std::vector<int> indices(vectorSize);
  std::iota (indices.begin(), indices.end(), 0);    // fill with increasing numbers: 0,1,2,...
  std::vector<double> values(vectorSize);
  VecGetValues(vector, vectorSize, indices.data(), values.data());

  return encodeBase64(values, withEncodedSizePrefix);
}

std::string Paraview::convertToAscii(const Vec &vector, bool fixedFormat)
{
  std::vector<double> vectorValues;
  PetscUtility::getVectorEntries(vector, vectorValues);

  return convertToAscii(vectorValues, fixedFormat);
}

std::string Paraview::encodeBase64(const std::vector<double> &vector, bool withEncodedSizePrefix)
{
  // encode as Paraview Float32
  assert(sizeof(float) == 4);

  int rawLength = vector.size()*sizeof(float);
  int dataStartPos = 0;

  if (withEncodedSizePrefix)
  {
    rawLength += 4;
    dataStartPos = 4;
  }

  int encodedLength = Base64::EncodedLength(rawLength);

  char raw[rawLength];
  // loop over vector entries and add bytes to raw buffer
  for (unsigned int i = 0; i < vector.size(); i++)
  {
    union {
      float d;    // note, this is a float, not a double, to get 4 bytes
      char c[4];
    };
    d = vector[i];
    memcpy(raw + dataStartPos + i*sizeof(float), c, 4);
  }

  // prepend number of bytes as uint32
  if (withEncodedSizePrefix)
  {
    union {
      uint32_t i;
      char c[4];
    };
    i = rawLength;

    memcpy(raw, c, 4);
  }

  char encoded[encodedLength+1];

  bool success = Base64::Encode(raw, rawLength, encoded, encodedLength);
  if (!success)
    LOG(WARNING) << "Base64 encoding failed";

  encoded[encodedLength] = '\0';

  return std::string(encoded);
}

std::string Paraview::encodeBase64(const std::vector<int> &vector, bool withEncodedSizePrefix)
{
  // encode as Paraview Int32
  assert(sizeof(int) == 4);

  int rawLength = vector.size()*sizeof(int);
  int dataStartPos = 0;

  if (withEncodedSizePrefix)
  {
    rawLength += 4;
    dataStartPos = 4;
  }

  int encodedLength = Base64::EncodedLength(rawLength);

  char raw[rawLength];
  // loop over vector entries and add bytes to raw buffer
  for (unsigned int i = 0; i < vector.size(); i++)
  {
    union {
      int integer;
      char c[4];
    };
    integer = vector[i];
    memcpy(raw + dataStartPos + i*sizeof(float), c, 4);
  }

  // prepend number of bytes as uint32
  if (withEncodedSizePrefix)
  {
    union {
      uint32_t i;
      char c[4];
    };
    i = rawLength;

    memcpy(raw, c, 4);
  }

  char encoded[encodedLength+1];

  bool success = Base64::Encode(raw, rawLength, encoded, encodedLength);
  if (!success)
    LOG(WARNING) << "Base64 encoding failed";


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
