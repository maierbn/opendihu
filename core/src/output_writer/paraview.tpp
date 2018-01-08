#include "output_writer/paraview.h"

#include <fstream>
#include <iostream>
#include <thread>
#include <chrono>

#include "easylogging++.h"
#include "base64.h"

#include <utility/python_utility.h>
#include <utility/petsc_utility.h>
#include <mesh/regular_fixed.h>
#include <mesh/structured_deformable.h>
#include <mesh/unstructured_deformable.h>
#include <mesh/mesh.h>

namespace OutputWriter
{

template<typename DataType>
void Paraview::writeSolution(DataType& data, int timeStepNo, double currentTime)
{
  if (!data.mesh())
  {
    LOG(FATAL) << "mesh is not set!";
  }
  const int dimension = data.mesh()->dimension();
  
  // solution and rhs vectors in mesh shape
  switch(dimension)
  {
  case 1:
    writeSolutionDim<1>(data, timeStepNo, currentTime);
    break;
  case 2:
    writeSolutionDim<2>(data, timeStepNo, currentTime);
    break;
  case 3:
    writeSolutionDim<3>(data, timeStepNo, currentTime);
    break;
  };
}

template <int dimension, typename DataType>
void Paraview::writeSolutionDim(DataType &data, int timeStepNo, double currentTime)
{
  LOG(TRACE) << "writeMesh<"<<dimension<<">()";
  
  if (std::dynamic_pointer_cast<Mesh::RegularFixed<dimension>>(data.mesh()) != NULL)
  {
    writeRectilinearGrid<Mesh::RegularFixed<dimension>>(data);
  }
  else if (std::dynamic_pointer_cast<
             BasisFunction::BasisOnMesh<Mesh::StructuredDeformable<dimension>,BasisFunction::Lagrange<1>>
           >(data.mesh()) != NULL)
  {
    // structured grid only for elements that contain only nodes at the corners (i.e. linear lagrange elements)
    writeStructuredGrid<dimension>(data);
  }
  else if (std::dynamic_pointer_cast<Mesh::UnstructuredDeformable<dimension>>(data.mesh()) != NULL)
  {
    writeUnstructuredGrid<dimension>(data);
  }
}

template <class Mesh>
void Paraview::writeRectilinearGrid(Data::Data& data)
{
  // determine file name
  std::stringstream s;
  s<<filename_<<".vtr";
  std::string filename = s.str();

  // open file
  std::ofstream file = openFile(filename);
  
  LOG(DEBUG) << "Write RectilinearGrid, file \""<<filename<<"\".";
  
  // determine values
  std::shared_ptr<Mesh> mesh = std::static_pointer_cast<Mesh>(data.mesh());
  
  // extent
  std::vector<int> extent = {0,0,0};   // number of nodes in x, y and z direction
  for (int i=0; i<mesh->dimension(); i++)
    extent[i] = mesh->nElements(i) + 1 - 1;   // nDOFS = nElements+1, value-1 because indices start with 0
  
  // coordinates of grid 
  std::vector<double> coordinates[3];
  int dimensionNo = 0;
  for (; dimensionNo<mesh->dimension(); dimensionNo++)
  {
    double meshWidth = mesh->meshWidth(dimensionNo);
    double nElements = mesh->nElements(dimensionNo);
    
    LOG(DEBUG) << "dimension "<<dimensionNo<<", meshWidth: "<<meshWidth<<", nElements: "<<nElements;
    
    coordinates[dimensionNo].resize(nElements+1);
    
    for(int nodeNo = 0; nodeNo < nElements+1; nodeNo++)
    {
      double coordinate = nodeNo * meshWidth;
      LOG(DEBUG) << "coordinate: "<<coordinate<<", nodeNo="<<nodeNo;
      coordinates[dimensionNo][nodeNo] = coordinate;
    }
  }
  
  // set other coordinates to 0
  for(; dimensionNo<3; dimensionNo++)
  {
    coordinates[dimensionNo].resize(1);
    coordinates[dimensionNo][0] = 0.0;
  }
  
  // name of value field
  std::string scalarsName = "Solution";
  
  bool binaryOutput = PythonUtility::getOptionBool(specificSettings_, "binaryOutput", true);
  bool fixedFormat = PythonUtility::getOptionBool(specificSettings_, "fixedFormat", true);
  
  // write file
  file << "<?xml version=\"1.0\"?>" << std::endl
    << "<VTKFile type=\"RectilinearGrid\" version=\"1.0\" byte_order=\"LittleEndian\">" << std::endl    // intel cpus are LittleEndian
    << std::string(1, '\t') << "<RectilinearGrid " 
      << "WholeExtent=\"" << "0 " << extent[0] << " 0 "<< extent[1] << " 0 " << extent[2] << "\"> " << std::endl     // dataset element
    << std::string(2, '\t') << "<Piece " 
      << "Extent=\"0 " << extent[0] << " 0 "<< extent[1] << " 0 " << extent[2] << "\"> " << std::endl
    << std::string(3, '\t') << "<PointData Scalars=\"" << scalarsName << "\">" << std::endl;
  if (binaryOutput)
  {
    file << std::string(4, '\t') << "<DataArray " 
        << "Name=\"" << scalarsName << "\" " 
        << "type=\"Float32\" " 
        << "NumberOfComponents=\"1\" "
        << "format=\"binary\" >" << std::endl
      << std::string(5, '\t') << encodeBase64(data.solution()) << std::endl
      << std::string(4, '\t') << "</DataArray>" << std::endl;
  }
  else
  {
    file << std::string(4, '\t') << "<DataArray " 
        << "Name=\"" << scalarsName << "\" " 
        << "type=\"Float32\" " 
        << "NumberOfComponents=\"1\" "
        << "format=\"ascii\" >" << std::endl
      << std::string(5, '\t') << convertToAscii(data.solution(), fixedFormat) << std::endl
      << std::string(4, '\t') << "</DataArray>" << std::endl;
  }
  
  file << std::string(3, '\t') << "</PointData>" << std::endl
    << std::string(3, '\t') << "<CellData>" << std::endl
    << std::string(3, '\t') << "</CellData>" << std::endl
    << std::string(3, '\t') << "<Coordinates>" << std::endl;
    
  if (binaryOutput)
  {
    file << std::string(4, '\t') << "<DataArray " 
        << "type=\"Float32\" " 
        << "NumberOfComponents=\"1\" "
        << "format=\"binary\" >" << std::endl
      << std::string(5, '\t') << encodeBase64(coordinates[0]) << std::endl
      << std::string(4, '\t') << "</DataArray>" << std::endl
      << std::string(4, '\t') << "<DataArray " 
        << "type=\"Float32\" " 
        << "NumberOfComponents=\"1\" "
        << "format=\"binary\" >" << std::endl
      << std::string(5, '\t') << encodeBase64(coordinates[1]) << std::endl
      << std::string(4, '\t') << "</DataArray>" << std::endl
      << std::string(4, '\t') << "<DataArray " 
        << "type=\"Float32\" " 
        << "NumberOfComponents=\"1\" "
        << "format=\"binary\" >" << std::endl
      << std::string(5, '\t') << encodeBase64(coordinates[2]) << std::endl
      << std::string(4, '\t') << "</DataArray>" << std::endl;
  }
  else
  {
    file << std::string(4, '\t') << "<DataArray " 
        << "type=\"Float32\" " 
        << "NumberOfComponents=\"1\" "
        << "format=\"ascii\" >" << std::endl
      << std::string(5, '\t') << convertToAscii(coordinates[0], fixedFormat) << std::endl
      << std::string(4, '\t') << "</DataArray>" << std::endl
      << std::string(4, '\t') << "<DataArray " 
        << "type=\"Float64\" " 
        << "NumberOfComponents=\"1\" "
        << "format=\"ascii\" >" << std::endl
      << std::string(5, '\t') << convertToAscii(coordinates[1], fixedFormat) << std::endl
      << std::string(4, '\t') << "</DataArray>" << std::endl
      << std::string(4, '\t') << "<DataArray " 
        << "type=\"Float64\" " 
        << "NumberOfComponents=\"1\" "
        << "format=\"ascii\" >" << std::endl
      << std::string(5, '\t') << convertToAscii(coordinates[2], fixedFormat) << std::endl
      << std::string(4, '\t') << "</DataArray>" << std::endl;
  }
  file << std::string(3, '\t') << "</Coordinates>" << std::endl
    << std::string(2, '\t') << "</Piece>" << std::endl
    << std::string(1, '\t') << "</RectilinearGrid>" << std::endl
    << "</VTKFile>"<<std::endl;
}

template <int D>
void Paraview::writeStructuredGrid(Data::Data& data)
{
  // determine file name
  std::stringstream s;
  s<<filename_<<".vts";
  std::string filename = s.str();

  // open file
  std::ofstream file = openFile(filename);
  
  LOG(DEBUG) << "Write StructuredGrid, file \""<<filename<<"\".";
  
  std::shared_ptr<Mesh::StructuredDeformable<D>> mesh = std::static_pointer_cast<Mesh::StructuredDeformable<D>>(data.mesh());
  
  // extent
  std::vector<int> extent = {0,0,0};   // number of nodes in x, y and z direction
  for (int i=0; i<mesh->dimension(); i++)
    extent[i] = mesh->nElements(i) + 1 - 1;   // nDOFS = nElements+1, value-1 because indices start with 0
  
  // points
  std::vector<double> points;
  mesh->getNodePositions(points);
  
  // name of value field
  std::string scalarsName = "Solution";
  
  bool binaryOutput = PythonUtility::getOptionBool(specificSettings_, "binaryOutput", true);
  bool fixedFormat = PythonUtility::getOptionBool(specificSettings_, "fixedFormat", true);
  
  // write file
  file << "<?xml version=\"1.0\"?>" << std::endl
    << "<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">" << std::endl    // intel cpus are LittleEndian
    << std::string(1, '\t') << "<StructuredGrid " 
      << "WholeExtent=\"" << "0 " << extent[0] << " 0 "<< extent[1] << " 0 " << extent[2] << "\"> " << std::endl     // dataset element
    << std::string(2, '\t') << "<Piece " 
      << "Extent=\"0 " << extent[0] << " 0 "<< extent[1] << " 0 " << extent[2] << "\"> " << std::endl
    << std::string(3, '\t') << "<PointData Scalars=\"" << scalarsName << "\">" << std::endl;
  if (binaryOutput)
  {
    file << std::string(4, '\t') << "<DataArray " 
        << "Name=\"" << scalarsName << "\" " 
        << "type=\"Float32\" " 
        << "NumberOfComponents=\"1\" "
        << "format=\"binary\" >" << std::endl
      << std::string(5, '\t') << encodeBase64(data.solution()) << std::endl
      << std::string(4, '\t') << "</DataArray>" << std::endl;
  }
  else
  {
    file << std::string(4, '\t') << "<DataArray " 
        << "Name=\"" << scalarsName << "\" " 
        << "type=\"Float32\" " 
        << "NumberOfComponents=\"1\" "
        << "format=\"ascii\" >" << std::endl
      << std::string(5, '\t') << convertToAscii(data.solution(), fixedFormat) << std::endl
      << std::string(4, '\t') << "</DataArray>" << std::endl;
  }
  
  file << std::string(3, '\t') << "</PointData>" << std::endl
    << std::string(3, '\t') << "<CellData>" << std::endl
    << std::string(3, '\t') << "</CellData>" << std::endl
    << std::string(3, '\t') << "<Points>" << std::endl;
    
  if (binaryOutput)
  {
    file << std::string(4, '\t') << "<DataArray " 
        << "type=\"Float32\" " 
        << "NumberOfComponents=\"3\" "
        << "format=\"binary\" >" << std::endl
      << std::string(5, '\t') << encodeBase64(points) << std::endl
      << std::string(4, '\t') << "</DataArray>" << std::endl;
  }
  else
  {
    file << std::string(4, '\t') << "<DataArray " 
        << "type=\"Float32\" " 
        << "NumberOfComponents=\"3\" "
        << "format=\"ascii\" >" << std::endl
      << std::string(5, '\t') << convertToAscii(points, fixedFormat) << std::endl
      << std::string(4, '\t') << "</DataArray>" << std::endl;
  }
  file << std::string(3, '\t') << "</Points>" << std::endl
    << std::string(2, '\t') << "</Piece>" << std::endl
    << std::string(1, '\t') << "</StructuredGrid>" << std::endl
    << "</VTKFile>"<<std::endl;
}

template <int D>
void Paraview::writeUnstructuredGrid(Data::Data& data)
{
  // determine file name
  std::stringstream s;
  s<<filename_<<".vtu";
  std::string filename = s.str();

  // open file
  std::ofstream file = openFile(filename);
  
  LOG(DEBUG) << "Write UnstructuredGrid, file \""<<filename<<"\".";
  
  std::shared_ptr<Mesh::Deformable<D>> mesh = std::static_pointer_cast<Mesh::Deformable<D>>(data.mesh());
  
  // extent
  std::vector<int> extent = {0,0,0};   // number of nodes in x, y and z direction
  for (int i=0; i<mesh->dimension(); i++)
    extent[i] = mesh->nElements(i) + 1 - 1;   // nDOFS = nElements+1, value-1 because indices start with 0
  
  // points
  std::vector<double> points;
  mesh->getNodePositions(points);
  
  // name of value field
  std::string scalarsName = "Solution";
  
  bool binaryOutput = PythonUtility::getOptionBool(specificSettings_, "binaryOutput", true);
  bool fixedFormat = PythonUtility::getOptionBool(specificSettings_, "fixedFormat", true);
  
  // write file
  file << "<?xml version=\"1.0\"?>" << std::endl
    << "<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">" << std::endl    // intel cpus are LittleEndian
    << std::string(1, '\t') << "<StructuredGrid " 
      << "WholeExtent=\"" << "0 " << extent[0] << " 0 "<< extent[1] << " 0 " << extent[2] << "\"> " << std::endl     // dataset element
    << std::string(2, '\t') << "<Piece " 
      << "Extent=\"0 " << extent[0] << " 0 "<< extent[1] << " 0 " << extent[2] << "\"> " << std::endl
    << std::string(3, '\t') << "<PointData Scalars=\"" << scalarsName << "\">" << std::endl;
  if (binaryOutput)
  {
    file << std::string(4, '\t') << "<DataArray " 
        << "Name=\"" << scalarsName << "\" " 
        << "type=\"Float32\" " 
        << "NumberOfComponents=\"1\" "
        << "format=\"binary\" >" << std::endl
      << std::string(5, '\t') << encodeBase64(data.solution()) << std::endl
      << std::string(4, '\t') << "</DataArray>" << std::endl;
  }
  else
  {
    file << std::string(4, '\t') << "<DataArray " 
        << "Name=\"" << scalarsName << "\" " 
        << "type=\"Float32\" " 
        << "NumberOfComponents=\"1\" "
        << "format=\"ascii\" >" << std::endl
      << std::string(5, '\t') << convertToAscii(data.solution(), fixedFormat) << std::endl
      << std::string(4, '\t') << "</DataArray>" << std::endl;
  }
  
  file << std::string(3, '\t') << "</PointData>" << std::endl
    << std::string(3, '\t') << "<CellData>" << std::endl
    << std::string(3, '\t') << "</CellData>" << std::endl
    << std::string(3, '\t') << "<Points>" << std::endl;
    
  if (binaryOutput)
  {
    file << std::string(4, '\t') << "<DataArray " 
        << "type=\"Float32\" " 
        << "NumberOfComponents=\"3\" "
        << "format=\"binary\" >" << std::endl
      << std::string(5, '\t') << encodeBase64(points) << std::endl
      << std::string(4, '\t') << "</DataArray>" << std::endl;
  }
  else
  {
    file << std::string(4, '\t') << "<DataArray " 
        << "type=\"Float32\" " 
        << "NumberOfComponents=\"3\" "
        << "format=\"ascii\" >" << std::endl
      << std::string(5, '\t') << convertToAscii(points, fixedFormat) << std::endl
      << std::string(4, '\t') << "</DataArray>" << std::endl;
  }
  file << std::string(3, '\t') << "</Points>" << std::endl
    << std::string(2, '\t') << "</Piece>" << std::endl
    << std::string(1, '\t') << "</StructuredGrid>" << std::endl
    << "</VTKFile>"<<std::endl;
}

};