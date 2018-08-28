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
  int encodedLength = Base64::EncodedLength(rawLength);

  char raw[rawLength];
  for (unsigned int i=0; i<vector.size(); i++)
  {
    union {
      float d;
      char c[4];
    };
    d = vector[i];
    memcpy(raw+i*sizeof(float), c, 4);
  }

  char encoded[encodedLength+1];
  //Base64::Encode(reinterpret_cast<char *>(vector.data()), rawLength, encoded, encodedLength);
  bool success = Base64::Encode(raw, rawLength, encoded, encodedLength);
  if (!success)
    LOG(WARNING) << "encoding failed";

  encoded[encodedLength] = '\0';

  return std::string(encoded);
}

std::string Paraview::encodeBase64(const std::vector<element_no_t> &vector)
{
  // encode as Paraview Int32
  assert(sizeof(element_no_t) == 4);
  
  int rawLength = vector.size()*sizeof(element_no_t);
  int encodedLength = Base64::EncodedLength(rawLength);

  char raw[rawLength];
  for (unsigned int i=0; i<vector.size(); i++)
  {
    union {
      element_no_t integer;
      char c[4];
    };
    integer = vector[i];
    memcpy(raw+i*sizeof(element_no_t), c, 4);
  }

  char encoded[encodedLength+1];
  //Base64::Encode(reinterpret_cast<char *>(vector.data()), rawLength, encoded, encodedLength);
  bool success = Base64::Encode(raw, rawLength, encoded, encodedLength);
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

std::string Paraview::convertToAscii(const std::vector<element_no_t> &vector, bool fixedFormat)
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
/*
void Paraview::writeVTKMasterFile()
{

  file<< "<?xml version=\"1.0\"?>" <<endl
    << "<VTKFile type=\"PRectilinearGrid\">" <<endl
    << "<PRectilinearGrid WholeExtent=\"0 " << comm->get_global_dimension()[0]<< " 0 "
      << comm->get_global_dimension()[1]<< " 0 0\" GhostLevel=\"0\">"
      <<endl
    << "<PCoordinates>" <<endl
    << "<PDataArray type=\"Float64\"/>" <<endl
    << "<PDataArray type=\"Float64\"/>" <<endl
    << "<PDataArray type=\"Float64\"/>" <<endl
    << "</PCoordinates>" <<endl;

    MultiIndexType n_subareas = comm->get_n_subareas();

    int x_begin = 0;
    int y_begin = 0;

    MultiIndexType local_dimension;

    for(int y = 0; y < n_subareas[1]; y++)      // iterate over row (y-direction)
    {
        x_begin = 0;
        for(int x = 0; x < n_subareas[0]; x++)      // iterate over column (x-direction)
        {
            local_dimension = comm->get_local_dimension(x, y);

            os<< "<Piece Extent=\""
              <<x_begin<< " " <<x_begin+local_dimension[0]<< " "
              <<y_begin<< " " <<y_begin+local_dimension[1]<< " 0 0\" "
              << "Source=\"field_" <<step<< "_processor_" <<y<< "_" <<x<< ".vtr\"/>" <<endl;

            x_begin += local_dimension[0];
        }
        y_begin += local_dimension[1];
    }

  os<< "<PPointData Vectors=\"field\" Scalars=\"p, vorticity, stream\">" <<endl;
  os<< "<PDataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"field\"     format=\"ascii\"/>" <<endl;
  os<< "<PDataArray type=\"Float64\"                          Name=\"p\"         format=\"ascii\"/>" <<endl;
  os<< "<PDataArray type=\"Float64\"                          Name=\"vorticity\" format=\"ascii\"/>" <<endl;
  os<< "<PDataArray type=\"Float64\"                          Name=\"stream\"    format=\"ascii\"/>" <<endl;
    os<< "<PDataArray type=\"Float64\"                          Name=\"rhs\"       format=\"ascii\"/>" <<endl;
    os<< "<PDataArray type=\"Float64\"                          Name=\"partitioning\"    format=\"ascii\"/>" <<endl;
  os<< "</PPointData>" <<endl;

  os<< "</PRectilinearGrid>" <<endl;

    os<< "<PUnstructuredGrid GhostLevel=\"0\">" <<endl;
    os<< "<PPointData></PPointData>" <<endl;
    os<< "<PCellData>" <<endl;
    os<< "<DataArray type=\"Float32\" Name=\"particles\" format=\"ascii\">" <<endl;
    os<< "</DataArray>" <<endl;
    os<< "</PCellData>" <<endl;
    os<< "<PPoints>" <<endl;
    os<< "<PDataArray NumberOfComponents=\"3\">" <<endl;
    os<< "</PDataArray>" <<endl;
    os<< "</PPoints>" <<endl;
    os<< "<Piece Source=\"particles_" <<step<< ".vtu\"/>" <<endl;

    os<< "</PUnstructuredGrid> " <<endl;

  os<< "</VTKFile>" <<endl;
}

void Paraview::writeVTKSlaveFile()
{

}
*/

};