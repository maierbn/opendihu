#include "output_writer/paraview/paraview_writer.h"

#include <fstream>
#include <iostream>
#include <thread>
#include <chrono>

#include "output_writer/paraview/loop_collect_field_variables_names.h"
#include "output_writer/paraview/loop_output_point_data.h"

namespace OutputWriter
{
// regular fixed
template<int D, typename BasisFunctionType, typename OutputFieldVariablesType>
void ParaviewWriter<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<D>, BasisFunctionType>, OutputFieldVariablesType>::
outputFile(std::string filename, OutputFieldVariablesType fieldVariables, std::string meshName, 
           std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<D>, BasisFunctionType>> mesh,
           int nFieldVariablesOfMesh, PythonConfig specificSettings)
{
  // write a RectilinearGrid

  // get type of geometry field
  typedef FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<D>, BasisFunctionType>,3> GeometryFieldType;

  // collect field variable names that are defined on the current mesh
  std::vector<std::string> namesScalars, namesVectors;
  ParaviewLoopOverTuple::loopCollectFieldVariablesNames(fieldVariables, meshName, namesScalars, namesVectors);

  if (specificSettings.hasKey("binaryOutput"))
  {
    LOG(ERROR) << "Key \"binaryOutput\" for Paraview output was recently changed to \"binary\"!";
  }
  bool binaryOutput = specificSettings.getOptionBool("binary", true);
  bool fixedFormat = specificSettings.getOptionBool("fixedFormat", true);

  // determine file name
  std::stringstream s;
  if (mesh->meshPartition()->nRanks() > 1 && mesh->meshPartition()->ownRankNo() == 0)
  {
    // this is a parallel run, write master file, "parallel dataset element"

    // extract filename base
    std::string filenameBase, filenameBaseWithPath;
    std::size_t pos1 = filename.rfind(".");
    if (pos1 != std::string::npos)
    {
      filenameBaseWithPath = filename.substr(0, pos1);
    }
    else
    {
      filenameBaseWithPath = filename;
    }

    pos1 = filename.rfind("/");
    if (pos1 != std::string::npos)
    {
      filenameBase = filenameBaseWithPath.substr(pos1+1);
    }
    else
    {
      filenameBase = filenameBaseWithPath;
    }

    s << filenameBaseWithPath << ".pvtr";

    // open file
    std::ofstream file = Paraview::openFile(s.str());

    LOG(DEBUG) << "Write PRectilinearGrid, file \"" << s.str() << "\".";

    std::array<node_no_t,6> globalExtent = {0};
    for (int dimensionNo = 0; dimensionNo < D; dimensionNo++)
    {
      globalExtent[dimensionNo] = mesh->meshPartition()->nNodesGlobal(dimensionNo) - 1;
    }

    // write file
    file << "<?xml version=\"1.0\"?>" << std::endl
      << "<!-- " << DihuContext::versionText() << " " << DihuContext::metaText() << "-->" << std::endl
      << "<VTKFile type=\"PRectilinearGrid\" version=\"1.0\" byte_order=\"LittleEndian\">" << std::endl    // intel cpus are LittleEndian
      << std::string(1, '\t') << "<PRectilinearGrid "
      << "WholeExtent=\"" << "0 " << globalExtent[0] << " 0 " << globalExtent[1] << " 0 " << globalExtent[2] << "\" GhostLevel=\"0\"> " << std::endl;

    file << std::string(2, '\t') << "<PPointData ";
    // output first name of scalar fields, this is the default field to be displayed
    if (!namesScalars.empty())
    {
      file << "Scalars=\"" << namesScalars[0] << "\" ";
    }
      // output first name of vector fields, this is the default field to be displayed
    if (!namesVectors.empty())
    {
      file << "Vectors=\"" << namesVectors[0] << "\" ";
    }
    file << ">" << std::endl;

    ParaviewLoopOverTuple::loopOutputPointData(fieldVariables, meshName, file, binaryOutput, fixedFormat, true);
    Paraview::writeParaviewPartitionFieldVariable<GeometryFieldType>(mesh->geometryField(), file, binaryOutput, fixedFormat, true);

    file << std::string(2, '\t') << "</PPointData>" << std::endl
      << std::string(2, '\t') << "<PCellData>" << std::endl
      << std::string(2, '\t') << "</PCellData>" << std::endl
      << std::string(2, '\t') << "<PCoordinates>" << std::endl;

    std::string format;
    if (binaryOutput)
    {
      format = "binary";
    }
    else
    {
      format = "ascii";
    }
    file << std::string(3, '\t') << "<PDataArray "
      << "type=\"Float32\" "
      << "NumberOfComponents=\"1\" "
      << "format=\"" << format << "\" />" << std::endl
      << std::string(3, '\t') << "<PDataArray "
      << "type=\"Float32\" "
      << "NumberOfComponents=\"1\" "
      << "format=\"" << format << "\" />" << std::endl
      << std::string(3, '\t') << "<PDataArray "
      << "type=\"Float32\" "
      << "NumberOfComponents=\"1\" "
      << "format=\"" << format << "\" />" << std::endl
      << std::string(2, '\t') << "</PCoordinates>" << std::endl;

    for (int rankNo = 0; rankNo < mesh->meshPartition()->nRanks(); rankNo++)
    {
      // extent is number nodes-1 (not number elements)
      std::array<node_no_t,6> extent = {0};
      for (int dimensionNo = 0; dimensionNo < D; dimensionNo++)
      {
        int partitionIndex = mesh->meshPartition()->convertRankNoToPartitionIndex(dimensionNo, rankNo);
        extent[2*dimensionNo + 0] = mesh->meshPartition()->beginNodeGlobalNatural(dimensionNo, partitionIndex);
        extent[2*dimensionNo + 1] = extent[2*dimensionNo + 0] + mesh->meshPartition()->nNodesLocalWithGhosts(dimensionNo, partitionIndex)-1;
      }

      LOG(DEBUG) << "extent: " << extent;

      file << std::string(2, '\t') << "<Piece Extent=\"" << extent[0];
      for (int i = 1; i < 6; i++)
        file << " " << extent[i];

      std::stringstream pieceFilename;
      Generic::appendRankNo(pieceFilename, mesh->meshPartition()->nRanks(), rankNo);
      file << "\" Source=\"" << filenameBase << pieceFilename.str() << ".vtr\" />" << std::endl;
    }

    file << std::string(1, '\t') << "</PRectilinearGrid>" << std::endl
      << "</VTKFile>" << std::endl;

    file.close();

    // write serial slave file
    s.str("");
    s << filename << ".vtr";
  }
  else
  {
    s << filename << ".vtr";
  }


  // open file
  std::ofstream file = Paraview::openFile(s.str());

  LOG(DEBUG) << "Write RectilinearGrid, file \"" << s.str() << "\".";

  // extent
  std::array<node_no_t,6> localExtent = {0};
  std::array<node_no_t,6> globalExtent = {0};
  for (int dimensionNo = 0; dimensionNo < D; dimensionNo++)
  {
    localExtent[2*dimensionNo + 0] = mesh->meshPartition()->beginNodeGlobalNatural(dimensionNo);
    localExtent[2*dimensionNo + 1] = localExtent[2*dimensionNo + 0] + mesh->meshPartition()->nNodesLocalWithGhosts(dimensionNo)-1;
    globalExtent[dimensionNo] = mesh->meshPartition()->nNodesGlobal(dimensionNo) - 1;
  }

  // coordinates of grid
  std::array<std::vector<double>,3> coordinates;
  int dimensionNo = 0;
  for (; dimensionNo < D; dimensionNo++)
  {
    double meshWidth = mesh->meshWidth();
    node_no_t nNodes = mesh->meshPartition()->nNodesLocalWithGhosts(dimensionNo);
    
    LOG(DEBUG) << "dimension " << dimensionNo << ", meshWidth: " << meshWidth;

    coordinates[dimensionNo].resize(nNodes);

    for (node_no_t nodeNo = 0; nodeNo < nNodes; nodeNo++)
    {
      double coordinate = (mesh->meshPartition()->beginNodeGlobalNatural(dimensionNo) + nodeNo) * meshWidth;
      VLOG(1) << "coordinate: " << coordinate << ", nodeNo=" << nodeNo;
      coordinates[dimensionNo][nodeNo] = coordinate;
    }
  }

  // set other coordinates to 0
  for (; dimensionNo < 3; dimensionNo++)
  {
    coordinates[dimensionNo].resize(1);
    coordinates[dimensionNo][0] = 0.0;
  }
  
  // write file
  file << "<?xml version=\"1.0\"?>" << std::endl
    << "<!-- " << DihuContext::versionText() << " " << DihuContext::metaText() << "-->" << std::endl
    << "<VTKFile type=\"RectilinearGrid\" version=\"1.0\" byte_order=\"LittleEndian\">" << std::endl    // intel cpus are LittleEndian
    << std::string(1, '\t') << "<RectilinearGrid "
    << "WholeExtent=\"" << "0 " << globalExtent[0] << " 0 " << globalExtent[1] << " 0 " << globalExtent[2] << "\"> " << std::endl     // dataset element
    << std::string(2, '\t') << "<Piece Extent=\"" << localExtent[0];

  for (int i = 1; i < 6; i++)
    file << " " << localExtent[i];

  file << "\"> " << std::endl;

  file << std::string(3, '\t') << "<PointData ";
  // output first name of scalar fields, this is the default field to be displayed
  if (!namesScalars.empty())
  {
    file << "Scalars=\"" << namesScalars[0] << "\" ";
  }
    // output first name of vector fields, this is the default field to be displayed
  if (!namesVectors.empty())
  {
    file << "Vectors=\"" << namesVectors[0] << "\" ";
  }
  file << ">" << std::endl;
    
  ParaviewLoopOverTuple::loopOutputPointData(fieldVariables, meshName, file, binaryOutput, fixedFormat, false);
  Paraview::writeParaviewPartitionFieldVariable<GeometryFieldType>(mesh->geometryField(), file, binaryOutput, fixedFormat, false);
  
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
      << std::string(5, '\t') << Paraview::encodeBase64Float(coordinates[0].begin(), coordinates[0].end()) << std::endl
      << std::string(4, '\t') << "</DataArray>" << std::endl
      << std::string(4, '\t') << "<DataArray "
        << "type=\"Float32\" "
        << "NumberOfComponents=\"1\" "
        << "format=\"binary\" >" << std::endl
      << std::string(5, '\t') << Paraview::encodeBase64Float(coordinates[1].begin(), coordinates[1].end()) << std::endl
      << std::string(4, '\t') << "</DataArray>" << std::endl
      << std::string(4, '\t') << "<DataArray "
        << "type=\"Float32\" "
        << "NumberOfComponents=\"1\" "
        << "format=\"binary\" >" << std::endl
      << std::string(5, '\t') << Paraview::encodeBase64Float(coordinates[2].begin(), coordinates[2].end()) << std::endl
      << std::string(4, '\t') << "</DataArray>" << std::endl;
  }
  else
  {
    file << std::string(4, '\t') << "<DataArray "
        << "type=\"Float32\" "
        << "NumberOfComponents=\"1\" "
        << "format=\"ascii\" >" << std::endl
      << std::string(5, '\t') << Paraview::convertToAscii(coordinates[0], fixedFormat) << std::endl
      << std::string(4, '\t') << "</DataArray>" << std::endl
      << std::string(4, '\t') << "<DataArray "
        << "type=\"Float32\" "
        << "NumberOfComponents=\"1\" "
        << "format=\"ascii\" >" << std::endl
      << std::string(5, '\t') << Paraview::convertToAscii(coordinates[1], fixedFormat) << std::endl
      << std::string(4, '\t') << "</DataArray>" << std::endl
      << std::string(4, '\t') << "<DataArray "
        << "type=\"Float32\" "
        << "NumberOfComponents=\"1\" "
        << "format=\"ascii\" >" << std::endl
      << std::string(5, '\t') << Paraview::convertToAscii(coordinates[2], fixedFormat) << std::endl
      << std::string(4, '\t') << "</DataArray>" << std::endl;
  }
  file << std::string(3, '\t') << "</Coordinates>" << std::endl
    << std::string(2, '\t') << "</Piece>" << std::endl
    << std::string(1, '\t') << "</RectilinearGrid>" << std::endl
    << "</VTKFile>" << std::endl;
  
}
  
// structured deformable
template<int D, typename BasisFunctionType, typename OutputFieldVariablesType>
void ParaviewWriter<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>, BasisFunctionType>, OutputFieldVariablesType>::
outputFile(std::string filename, OutputFieldVariablesType fieldVariables, std::string meshName, 
           std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>, BasisFunctionType>> mesh,
           int nFieldVariablesOfMesh, PythonConfig specificSettings)
{
  // write a StructuredGrid

  // get type of geometry field
  typedef FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>, BasisFunctionType>,3> GeometryFieldType;

  // collect field variable names that are defined on the current mesh
  std::vector<std::string> namesScalars, namesVectors;
  ParaviewLoopOverTuple::loopCollectFieldVariablesNames(fieldVariables, meshName, namesScalars, namesVectors);

  if (specificSettings.hasKey("binaryOutput"))
  {
    LOG(ERROR) << "Key \"binaryOutput\" for Paraview output was recently changed to \"binary\"!";
  }
  bool binaryOutput = specificSettings.getOptionBool("binary", true);
  bool fixedFormat = specificSettings.getOptionBool("fixedFormat", true);

  // determine file name
  std::stringstream s;
  if (mesh->meshPartition()->nRanks() > 1 && mesh->meshPartition()->ownRankNo() == 0)
  {
    // this is a parallel run, write master file, "parallel dataset element"

    // extract filename base
    std::string filenameBase, filenameBaseWithPath;
    std::size_t pos1 = filename.rfind(".");
    if (pos1 != std::string::npos)
    {
      filenameBaseWithPath = filename.substr(0, pos1);
    }
    else
    {
      filenameBaseWithPath = filename;
    }

    pos1 = filename.rfind("/");
    if (pos1 != std::string::npos)
    {
      filenameBase = filenameBaseWithPath.substr(pos1+1);
    }
    else
    {
      filenameBase = filenameBaseWithPath;
    }

    s << filenameBaseWithPath << ".pvts";

    // open file
    std::ofstream file = Paraview::openFile(s.str());

    LOG(DEBUG) << "Write PStructuredGrid, file \"" << s.str() << "\".";

    std::array<node_no_t,6> globalExtent = {0};
    for (int dimensionNo = 0; dimensionNo < D; dimensionNo++)
    {
      globalExtent[dimensionNo] = mesh->meshPartition()->nNodesGlobal(dimensionNo) - 1;
    }

    // write file
    file << "<?xml version=\"1.0\"?>" << std::endl
      << "<!-- " << DihuContext::versionText() << " " << DihuContext::metaText() << "-->" << std::endl
      << "<VTKFile type=\"PStructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">" << std::endl    // intel cpus are LittleEndian
      << std::string(1, '\t') << "<PStructuredGrid "
      << "WholeExtent=\"" << "0 " << globalExtent[0] << " 0 " << globalExtent[1] << " 0 " << globalExtent[2] << "\" GhostLevel=\"0\"> " << std::endl;

    file << std::string(2, '\t') << "<PPointData ";
    // output first name of scalar fields, this is the default field to be displayed
    if (!namesScalars.empty())
    {
      file << "Scalars=\"" << namesScalars[0] << "\" ";
    }
      // output first name of vector fields, this is the default field to be displayed
    if (!namesVectors.empty())
    {
      file << "Vectors=\"" << namesVectors[0] << "\" ";
    }
    file << ">" << std::endl;

    ParaviewLoopOverTuple::loopOutputPointData(fieldVariables, meshName, file, binaryOutput, fixedFormat, true);
    Paraview::writeParaviewPartitionFieldVariable<GeometryFieldType>(mesh->geometryField(), file, binaryOutput, fixedFormat, true);

    file << std::string(2, '\t') << "</PPointData>" << std::endl
      << std::string(2, '\t') << "<PCellData>" << std::endl
      << std::string(2, '\t') << "</PCellData>" << std::endl
      << std::string(2, '\t') << "<PPoints>" << std::endl;

    std::string format;
    if (binaryOutput)
    {
      format = "binary";
    }
    else
    {
      format = "ascii";
    }
    file << std::string(3, '\t') << "<PDataArray "
      << "type=\"Float32\" "
      << "NumberOfComponents=\"3\" "
      << "format=\"" << format << "\" />" << std::endl
      << std::string(2, '\t') << "</PPoints>" << std::endl;

    for (int rankNo = 0; rankNo < mesh->meshPartition()->nRanks(); rankNo++)
    {
      // extent is number nodes-1 (not number elements)
      std::array<node_no_t,6> extent = {0};
      for (int dimensionNo = 0; dimensionNo < D; dimensionNo++)
      {
        int partitionIndex = mesh->meshPartition()->convertRankNoToPartitionIndex(dimensionNo, rankNo);
        extent[2*dimensionNo + 0] = mesh->meshPartition()->beginNodeGlobalNatural(dimensionNo, partitionIndex);
        extent[2*dimensionNo + 1] = extent[2*dimensionNo + 0] + mesh->meshPartition()->nNodesLocalWithGhosts(dimensionNo, partitionIndex)-1;
      }

      file << std::string(2, '\t') << "<Piece Extent=\"" << extent[0];
      for (int i = 1; i < 6; i++)
        file << " " << extent[i];

      std::stringstream pieceFilename;
      Generic::appendRankNo(pieceFilename, mesh->meshPartition()->nRanks(), rankNo);
      file << "\" Source=\"" << filenameBase << pieceFilename.str() << ".vts\" />" << std::endl;
    }

    file << std::string(1, '\t') << "</PStructuredGrid>" << std::endl
      << "</VTKFile>" << std::endl;

    file.close();

    // write serial slave file
    s.str("");
    s << filename << ".vts";
  }
  else
  {
    s << filename << ".vts";
  }

  // open file
  std::ofstream file = Paraview::openFile(s.str());

  LOG(DEBUG) << "Write StructuredGrid, file \"" << s.str() << "\".";

  // extent
  std::array<node_no_t,6> localExtent = {0};
  std::array<node_no_t,6> globalExtent = {0};
  for (int dimensionNo = 0; dimensionNo < D; dimensionNo++)
  {
    localExtent[2*dimensionNo + 0] = mesh->meshPartition()->beginNodeGlobalNatural(dimensionNo);
    localExtent[2*dimensionNo + 1] = localExtent[2*dimensionNo + 0] + mesh->meshPartition()->nNodesLocalWithGhosts(dimensionNo)-1;
    globalExtent[dimensionNo] = mesh->meshPartition()->nNodesGlobal(dimensionNo) - 1;
  }

  // avoid bug in paraview when reading binary encoded (base64) values for 1D meshes
  //if (D == 1 && extent[0] > 1 && binaryOutput)
  //   extent[0] -= 1;
  
  // write file
  file << "<?xml version=\"1.0\"?>" << std::endl
    << "<!-- " << DihuContext::versionText() << " " << DihuContext::metaText() << "-->" << std::endl
    << "<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">" << std::endl    // intel cpus are LittleEndian
    << std::string(1, '\t') << "<StructuredGrid "
    << "WholeExtent=\"" << "0 " << globalExtent[0] << " 0 " << globalExtent[1] << " 0 " << globalExtent[2] << "\"> " << std::endl     // dataset element
    << std::string(2, '\t') << "<Piece Extent=\"" << localExtent[0];

  for (int i = 1; i < 6; i++)
    file << " " << localExtent[i];

  file << "\"> " << std::endl;
  
  
  file << std::string(3, '\t') << "<PointData ";
  // output first name of scalar fields, this is the default field to be displayed
  if (!namesScalars.empty())
  {
    file << "Scalars=\"" << namesScalars[0] << "\" ";
  }

  // output first name of vector fields, this is the default field to be displayed
  if (!namesVectors.empty())
  {
    file << "Vectors=\"" << namesVectors[0] << "\" ";
  }
  file << ">" << std::endl;

  ParaviewLoopOverTuple::loopOutputPointData(fieldVariables, meshName, file, binaryOutput, fixedFormat, false);
  Paraview::writeParaviewPartitionFieldVariable<GeometryFieldType>(mesh->geometryField(), file, binaryOutput, fixedFormat, false);

  file << std::string(3, '\t') << "</PointData>" << std::endl
    << std::string(3, '\t') << "<CellData>" << std::endl
    << std::string(3, '\t') << "</CellData>" << std::endl
    << std::string(3, '\t') << "<Points>" << std::endl;

  Paraview::writeParaviewFieldVariable<GeometryFieldType>(mesh->geometryField(), file, binaryOutput, fixedFormat, false);

  file << std::string(3, '\t') << "</Points>" << std::endl
    << std::string(2, '\t') << "</Piece>" << std::endl
    << std::string(1, '\t') << "</StructuredGrid>" << std::endl
    << "</VTKFile>" << std::endl;
}
  
  
// unstructured deformable
template<int D, typename BasisFunctionType, typename OutputFieldVariablesType>
void ParaviewWriter<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, OutputFieldVariablesType>::
outputFile(std::string filename, OutputFieldVariablesType fieldVariables, std::string meshName, 
           std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>> mesh,
           int nFieldVariablesOfMesh, PythonConfig specificSettings)
{
  // write an UnstructuredGrid
  // determine file name
  std::stringstream s;
  s << filename << ".vtu";

  // open file
  std::ofstream file = Paraview::openFile(s.str());

  LOG(DEBUG) << "Write UnstructuredGrid, file \"" << s.str() << "\".";

  // get type of geometry field
  typedef FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType> FunctionSpace;
  typedef FieldVariable::FieldVariable<FunctionSpace,3> GeometryFieldType;
  
  
  // name of value field
  if (specificSettings.hasKey("binaryOutput"))
  {
    LOG(ERROR) << "Key \"binaryOutput\" for Paraview output was recently changed to \"binary\"!";
  }
  bool binaryOutput = specificSettings.getOptionBool("binary", true);
  bool fixedFormat = specificSettings.getOptionBool("fixedFormat", true);

  // write file
  file << "<?xml version=\"1.0\"?>" << std::endl
    << "<!-- " << DihuContext::versionText() << " " << DihuContext::metaText() << "-->" << std::endl
    << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">" << std::endl    // intel cpus are LittleEndian
    << std::string(1, '\t') << "<UnstructuredGrid> " << std::endl
    << std::string(2, '\t') << "<Piece "
    << "NumberOfPoints=\"" << mesh->nNodesGlobal() << "\" NumberOfCells=\"" << mesh->nElementsLocal() << "\">" << std::endl;
    
  // collect field variable names that are defined on the current mesh
  std::vector<std::string> namesScalars, namesVectors;
  ParaviewLoopOverTuple::loopCollectFieldVariablesNames(fieldVariables, meshName, namesScalars, namesVectors);
  
  
  file << std::string(3, '\t') << "<PointData ";
  // output first name of scalar fields, this is the default field to be displayed
  if (!namesScalars.empty())
  {
    file << "Scalars=\"" << namesScalars[0] << "\" ";
  }
    // output first name of vector fields, this is the default field to be displayed
  if (!namesVectors.empty())
  {
    file << "Vectors=\"" << namesVectors[0] << "\" ";
  }
  file << ">" << std::endl;
    
  ParaviewLoopOverTuple::loopOutputPointData(fieldVariables, meshName, file, binaryOutput, fixedFormat, false);
  Paraview::writeParaviewPartitionFieldVariable<GeometryFieldType>(mesh->geometryField(), file, binaryOutput, fixedFormat, false);

  file << std::string(3, '\t') << "</PointData>" << std::endl
    << std::string(3, '\t') << "<CellData>" << std::endl
    << std::string(3, '\t') << "</CellData>" << std::endl
    << std::string(3, '\t') << "<Points>" << std::endl;

  Paraview::writeParaviewFieldVariable<GeometryFieldType>(mesh->geometryField(), file, binaryOutput, fixedFormat);

  file << std::string(3, '\t') << "</Points>" << std::endl
    << std::string(3, '\t') << "<Cells>" << std::endl
    << std::string(4, '\t') << "<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" ";
    
  // get the elements point lists
  std::vector<node_no_t> values;
  values.reserve(mesh->nElementsLocal() * FunctionSpace::averageNNodesPerElement());
  
  // loop over elements and collect point numbers of the element
  for (element_no_t elementNo = 0; elementNo < mesh->nElementsLocal(); elementNo++)
  {
    std::array<dof_no_t,FunctionSpace::nDofsPerElement()> dofsOfElement = mesh->getElementDofNosLocal(elementNo);
    for (typename std::array<dof_no_t,FunctionSpace::nDofsPerElement()>::const_iterator iter = dofsOfElement.begin(); iter != dofsOfElement.end(); iter++)
    {
      dof_no_t dofNo = *iter;
      if (dofNo % FunctionSpace::nDofsPerNode() == 0)
      {
        node_no_t nodeNo = dofNo / FunctionSpace::nDofsPerNode();
        values.push_back(nodeNo);
      }
    }
  }  

  // write to file
  if (binaryOutput)
  {
    file << "format=\"binary\">" << std::endl
       << Paraview::encodeBase64Int32(values.begin(), values.end()) << std::endl;
  }
  else 
  {
    file << "format=\"ascii\">" << std::endl << std::string(5, '\t')
       << Paraview::convertToAscii(values, fixedFormat) << std::endl;
  }
  
  file << std::string(4, '\t') << "</DataArray>" << std::endl 
    << std::string(4, '\t') << "<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" ";
    
  // offsets 
  values.clear();
  values.resize(mesh->nElementsLocal());
  for (element_no_t elementNo = 0; elementNo < mesh->nElementsLocal(); elementNo++)
  {
    values[elementNo] = (elementNo + 1) * FunctionSpace::nNodesPerElement();
  }
    
  if (binaryOutput)
  {
    file << "format=\"binary\">" << std::endl
      << Paraview::encodeBase64Int32(values.begin(), values.end()) << std::endl;
  }
  else
  {
    file << "format=\"ascii\">" << std::endl << std::string(5, '\t')
      << Paraview::convertToAscii(values, fixedFormat) << std::endl;
  }
  
  file << std::string(4, '\t') << "</DataArray>" << std::endl 
    << std::string(4, '\t') << "<DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\">" << std::endl
    << std::string(5, '\t');
    
  // cell types
  int cellType = 0;
  switch (D)
  {
  case 1:
    cellType = 3; // VTK_LINE
    break;
  case 2:
    cellType = 8; // VTK_PIXEL
    break;
  case 3:
    cellType = 11; // VTK_VOXEL
    break;
  }
  for (element_no_t elementNo = 0; elementNo < mesh->nElementsLocal(); elementNo++)
  {
    file << cellType << " ";
  }
  file << std::endl << std::string(4, '\t') << "</DataArray>" << std::endl;
    
  file << std::string(3, '\t') << "</Cells>" << std::endl
    << std::string(2, '\t') << "</Piece>" << std::endl
    << std::string(1, '\t') << "</UnstructuredGrid>" << std::endl
    << "</VTKFile>" << std::endl;
}
  
}  // namespace
