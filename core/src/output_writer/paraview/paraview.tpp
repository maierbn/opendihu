#include "output_writer/paraview/paraview.h"

#include <fstream>
#include <iostream>
#include <thread>
#include <chrono>

#include "easylogging++.h"
#include "base64.h"

#include "output_writer/paraview/loop_output.h"
#include "output_writer/paraview/loop_collect_mesh_properties.h"
#include "output_writer/paraview/poly_data_properties_for_mesh.h"

namespace OutputWriter
{

template<typename DataType>
void Paraview::write(DataType& data, int timeStepNo, double currentTime)
{
  // check if output should be written in this timestep and prepare filename
  if (!Generic::prepareWrite(data, timeStepNo, currentTime))
  {
    return;
  }

  if (combineFile_)
  {
    Paraview::writePolyDataFile<typename DataType::OutputFieldVariables>(data.getOutputFieldVariables());
  }
  else
  {
    // collect all available meshes
    std::set<std::string> meshNames;
    LoopOverTuple::loopCollectMeshNames<typename DataType::OutputFieldVariables>(data.getOutputFieldVariables(), meshNames);

    // loop over meshes and create a paraview file for each
    for (std::string meshName : meshNames)
    {
      // setup name of file
      std::stringstream filenameStart;
      if (meshNames.size() == 1)
        filenameStart << this->filename_;
      else
        filenameStart << this->filename_ << "_" << meshName;

      // loop over all field variables and output those that are associated with the mesh given by meshName
      ParaviewLoopOverTuple::loopOutput(data.getOutputFieldVariables(), data.getOutputFieldVariables(), meshName, filenameStart.str(), specificSettings_);
    }
  }
}

template<typename FieldVariableType>
void Paraview::writeParaviewFieldVariable(FieldVariableType &fieldVariable, 
                                          std::ofstream &file, bool binaryOutput, bool fixedFormat, bool onlyParallelDatasetElement)
{
  // here we have the type of the mesh with meshName (which is typedef to FunctionSpace)
  //typedef typename FieldVariableType::FunctionSpace FunctionSpace;

  // if only the "parallel dataset element" stub which is needed in the master files, should be written
  if (onlyParallelDatasetElement)
  {
    file << std::string(3, '\t') << "<PDataArray "
        << "Name=\"" << fieldVariable.name() << "\" "
        << "type=\"Float32\" "
        << "NumberOfComponents=\"" << fieldVariable.nComponents() << "\" ";

    if (binaryOutput)
    {
      file << "format=\"binary\" />" << std::endl;
    }
    else
    {
      file << "format=\"ascii\" />" << std::endl;
    }
  }
  else
  {
    // write normal data element
    file << std::string(4, '\t') << "<DataArray "
        << "Name=\"" << fieldVariable.name() << "\" "
        << "type=\"Float32\" "
        << "NumberOfComponents=\"" << fieldVariable.nComponents() << "\" ";

    const int nComponents = FieldVariableType::nComponents();
    std::string stringData;

    std::vector<double> values;
    std::array<std::vector<double>, nComponents> componentValues;

    // initialize the dofNosLocalNaturalOrdering vector of the meshPartition to be able to get the values in the natural ordering
    fieldVariable.functionSpace()->meshPartition()->initializeDofNosLocalNaturalOrdering(fieldVariable.functionSpace());

    // ensure that ghost values are in place
    fieldVariable.startGhostManipulation();

    // get all local values including ghosts for the components
    for (int componentNo = 0; componentNo < nComponents; componentNo++)
    {
      std::vector<double> retrievedLocalValues;
      fieldVariable.getValues(componentNo, fieldVariable.functionSpace()->meshPartition()->dofNosLocalNaturalOrdering(),
                              retrievedLocalValues);

      const int nDofsPerNode = FieldVariableType::FunctionSpace::nDofsPerNode();
      const node_no_t nNodesLocal = fieldVariable.functionSpace()->meshPartition()->nNodesLocalWithGhosts();

      // for Hermite only extract the non-derivative values
      componentValues[componentNo].resize(nNodesLocal);

      int index = 0;
      for (int i = 0; i < nNodesLocal; i++)
      {
        componentValues[componentNo][i] = retrievedLocalValues[index];
        index += nDofsPerNode;
      }
    }
    values.reserve(componentValues[0].size()*nComponents);

    // copy values in consecutive order (x y z x y z) to output
    for (int i = 0; i < componentValues[0].size(); i++)
    {
      for (int componentNo = 0; componentNo < nComponents; componentNo++)
      {
        values.push_back(componentValues[componentNo][i]);
      }
    }

    if (binaryOutput)
    {
      stringData = Paraview::encodeBase64(values);
      file << "format=\"binary\" >" << std::endl;
    }
    else
    {
      stringData = Paraview::convertToAscii(values, fixedFormat);
      file << "format=\"ascii\" >" << std::endl;
    }
    file << std::string(5, '\t') << stringData << std::endl
      << std::string(4, '\t') << "</DataArray>" << std::endl;
  }
}
  
template<typename FieldVariableType>
void Paraview::writeParaviewPartitionFieldVariable(FieldVariableType &geometryField,
                                                   std::ofstream &file, bool binaryOutput, bool fixedFormat, bool onlyParallelDatasetElement)
{
  // if only the "parallel dataset element" stub which is needed in the master files, should be written
  if (onlyParallelDatasetElement)
  {
    file << std::string(3, '\t') << "<PDataArray "
        << "Name=\"partitioning\" "
        << "type=\"Int32\" "
        << "NumberOfComponents=\"1\" ";

    if (binaryOutput)
    {
      file << "format=\"binary\" />" << std::endl;
    }
    else
    {
      file << "format=\"ascii\" />" << std::endl;
    }
  }
  else
  {
    // write normal data element
    file << std::string(4, '\t') << "<DataArray "
        << "Name=\"partitioning\" "
        << "type=\"Int32\" "
        << "NumberOfComponents=\"1\" ";

    std::string stringData;

    // get own rank no
    int ownRankNoCommWorld;
    MPIUtility::handleReturnValue(MPI_Comm_rank(MPI_COMM_WORLD, &ownRankNoCommWorld));

    const node_no_t nNodesLocal = geometryField.functionSpace()->meshPartition()->nNodesLocalWithGhosts();

    std::vector<int> values(nNodesLocal, ownRankNoCommWorld);

    if (binaryOutput)
    {
      stringData = Paraview::encodeBase64(values);
      file << "format=\"binary\" >" << std::endl;
    }
    else
    {
      stringData = Paraview::convertToAscii(values, fixedFormat);
      file << "format=\"ascii\" >" << std::endl;
    }
    file << std::string(5, '\t') << stringData << std::endl
      << std::string(4, '\t') << "</DataArray>" << std::endl;
  }
}

template<typename OutputFieldVariablesType>
void Paraview::writePolyDataFile(const OutputFieldVariablesType &fieldVariables)
{
  // output a *.vtp file

  // collect the size data that is needed to compute offsets for parallel file output
  std::map<std::string, PolyDataPropertiesForMesh> meshProperties;
  ParaviewLoopOverTuple::loopCollectMeshProperties<OutputFieldVariablesType>(fieldVariables, meshProperties);

  /*
   PolyDataPropertiesForMesh:
  int dimensionality;    ///< D=1: object is a VTK "Line", D=2, D=3: object is a VTK "Poly"
  global_no_t nPoints;   ///< the number of points needed for representing the mesh
  global_no_t nCells;    ///< the number of VTK "cells", i.e. "Lines" or "Polys", which is the opendihu number of "elements"

  std::vector<std::pair<std::string,int>> pointDataArrays;   ///< <name,nComponents> of PointData DataArray elements
  */

  /* one VTKPiece is the XML element that will be output as <Piece></Piece>. It is created from one or multiple opendihu meshes
   */
  struct VTKPiece
  {
    std::vector<std::string> meshNamesCombinedMeshes;   ///< the meshNames of the combined meshes, or only one meshName if it is not a merged mesh
    PolyDataPropertiesForMesh properties;   ///< the properties of the merged mesh

    int nVTKLines = 0;      ///< number of Line elements in the Piece
    int nVTKPolys = 0;      ///< number of Poly elements in th Piece
    std::string firstScalarName;   ///< name of the first scalar field variable of the mesh
    std::string firstVectorName;   ///< name of the first non-scalar field variable of the mesh

    //! constructor, initialize nPoints and nCells to 0
    VTKPiece()
    {
      properties.nPoints = 0;
      properties.nCells = 0;
    }

    //! assign the correct values to nVTKLines and nVTKPolys if properties has been set
    void setVTKValues()
    {
      if (properties.dimensionality == 1)
        nVTKLines = properties.nCells;
      else
        nVTKPolys = properties.nCells;

      // set values for firstScalarName and firstVectorName from the values in pointDataArrays
      for (auto pointDataArray : pointDataArrays)
      {
        if (firstScalarName == "" && pointDataArray.second == 1)
          firstScalarName = pointDataArray.first;
        if (firstVectorName == "" && pointDataArray.second != 1)
          firstVectorName = pointDataArray.first;
      }
    }
  };

  // vtkPieces contains at the beginning all merged meshes (can be 0), as given by the "mergeMeshes" option. Afterwards there are the normal, not merged meshes.
  std::vector<VTKPiece> vtkPieces(mergeMeshes_.size());

  // parse the collected properties of the meshes that will be output to the file
  for (std::map<std::string,PolyDataPropertiesForMesh>::iterator iter = meshProperties.begin(); iter != meshProperties.end(); iter++)
  {
    std::string meshName = iter->first;
    const PolyDataPropertiesForMesh &properties = iter->second;

    // check if this mesh should be combined with other meshes
    bool mergeMesh = false;
    int i = 0;
    for (std::vector<std::vector<std::string>>::iterator mergeMeshesListIter = mergeMeshes_.begin(); mergeMeshesListIter != mergeMeshes_.end(); mergeMeshesListIter++, i++)
    {
      // check if mesh is contained in current list of merged meshes
      for (std::vector<std::string>::iterator mergeMeshesNamesIter = mergeMeshesListIter->begin(); mergeMeshesNamesIter != mergeMeshesListIter->end(); mergeMeshesNamesIter++)
      {
        if (*mergeMeshesNamesIter == meshName)
        {
          vtkPieces[i].meshNamesCombinedMeshes.append(meshName);
          mergeMesh = true;

          // check if meshes can be merged
          if (!vtkPieces[i].properties.pointDataArrays.empty())   // if properties are already assigned by an earlier mesh
          {
            if (vtkPieces[i].properties.dimensionality != iter->dimensionality)
            {
              LOG(ERROR) << "Meshes " << *mergeMeshesNamesIter << " cannot be merged. Dimensionality mismatches for "
                << meshName << " (is " << iter->dimensionality << " instead of " << vtkPieces[i].properties.dimensionality << ")";
              mergeMesh = false;
            }
            if (vtkPieces[i].properties.pointDataArrays.size() != iter->pointDataArrays.size())
            {
              LOG(ERROR) << "Meshes " << *mergeMeshesNamesIter << " cannot be merged. Number of field variables mismatches for "
                << meshName << " (is " << iter->pointDataArrays.size() << " instead of " << vtkPieces[i].properties.pointDataArrays.size() << ")";
              mergeMesh = false;
            }
            else
            {
              for (int j = 0; j < iter->pointDataArrays.size(); j++)
              {
                if (vtkPieces[i].properties.pointDataArrays[j]->first != iter->pointDataArrays[j]->first)  // if the name of the jth field variable is different
                {
                  LOG(ERROR) << "Meshes " << *mergeMeshesNamesIter << " cannot be merged. Field variable names mismatch for "
                    << meshName << " (there is \"" << vtkPieces[i].properties.pointDataArrays[j]->first << "\" instead of \"" << iter->pointDataArrays[j]->first << "\")";
                  mergeMesh = false;
                }
              }
            }

            if (mergeMesh)
            {
              LOG(DEBUG) << "Merge mesh " << meshName << " into " << vtkPieces[i].meshNamesCombinedMeshes << ", add " << iter->nPoints << " points, " << iter->nCells << " elements to "
                << vtkPieces[i].properties.nPoints << " points, " << vtkPieces[i].properties.nCells << " elements";

              vtkPieces[i].properties.nPoints += iter->nPoints;
              vtkPieces[i].properties.nCells += iter->nCells;
              vtkPieces[i].setVTKValues();
              break;
            }
          }
          else
          {
            // properties are not yet assigned
            vtkPieces[i].properties = *iter; // store properties
            vtkPieces[i].setVTKValues();
            break;
          }

        }
      }
    }

    // if mesh should be merged with others, the corresponding vtkPiece was already update in the previous loops
    if (mergeMesh)
    {
      continue;
    }

    // add vtkPieces for this mesh
    vtkPieces.emplace_back();
    vtkPieces.back().meshNamesCombinedMeshes.push_back(meshName);
    vtkPieces.back().properties = *iter;
    vtkPieces.back().setVTKValues();

    LOG(DEBUG) << "Add mesh \"" << meshName << "\" as vtk piece with " << vtkPieces.back().properties.nPoints << " points and " << vtkPieces.back().properties.nCells
      << " cells, D=" << vtkPieces.back().properties.dimensionality << ", nVTKLines: " << vtkPieces.back().nVTKLines  << ", nVTKPolys: " << vtkPieces.back().nVTKPolys;
  }

  std::shared_ptr<Partition::RankSubset> rankSubsetAllComputedInstances = this->context_.partitionManager()->rankSubsetForCollectiveOperations();

  int nOutputFileParts = 1;

  // loop over Pieces
  for (std::vector<VTKPiece>::iterator vtkPieceIter = vtkPieces.begin(); vtkPieceIter != vtkPieces.end(); vtkPieceIter++, outputFilePartNo++)
  {
    nOutputFileParts += 2 + vtkPieceIter->properties.pointDataArrays.size()
      + (vtkPieceIter->nVTKLines != 0? 1 : 0)
      + (vtkPieceIter->nVTKPolys != 0? 1 : 0);
  }

  // create the basic structure of the output file
  std::stringstream outputFileParts[nOutputFileParts];
  int outputFilePartNo = 0;
  outputFileParts[outputFilePartNo] << "<?xml version=\"1.0\"?>" << std::endl
    << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\">" << std::endl    // intel cpus are LittleEndian
    << std::string(1, '\t') << "<PolyData>" << std::endl;
  outputFilePartNo++;

  // loop over Pieces
  for (std::vector<VTKPiece>::iterator vtkPieceIter = vtkPieces.begin(); vtkPieceIter != vtkPieces.end(); vtkPieceIter++, outputFilePartNo++)
  {
    outputFileParts[outputFilePartNo] << std::string(2, '\t') << "<Piece NumberOfPoints=\"" << vtkPieceIter->properties.nPoints << "\" NumberOfVerts=\"0\" "
      << "NumberOfLines=\"" << vtkPieceIter->nVTKLines << "\" NumberOfStrips=\"0\" NumberOfPolys=\"" << vtkPieceIter->nVTKPolys << "\">" << std::endl
      << std::string(3, '\t') << "<PointData";
    if (vtkPieceIter->firstScalarName != "")
    {
      outputFileParts[outputFilePartNo] << " Scalars=\"" << vtkPieceIter->firstScalarName << "\";
    }
    if (vtkPieceIter->firstVectorName != "")
    {
      outputFileParts[outputFilePartNo] << " Vectors=\"" << vtkPieceIter->firstVectorName; << "\"";
    }
    outputFileParts[outputFilePartNo] << ">";

    // loop over field variables (PointData)
    for (std::vector<std::pair<std::string,int>>::iterator pointDataArrayIter = vtkPieceIter->properties.pointDataArrays.begin(); pointDataArrayIter != vtkPieceIter->properties.pointDataArrays.end(); pointDataArrayIter++)
    {
      // write normal data element
      outputFileParts[outputFilePartNo] << std::string(4, '\t') << "<DataArray "
          << "Name=\"" << pointDataArrayIter->first << "\" "
          << "type=\"Float32\" "
          << "NumberOfComponents=\"" << pointDataArrayIter->second << "\" format=\"" << (binaryOutput_? "binary" : "ascii")
          << "\" >" << std::endl << std::string(5, '\t');

      // at this point the data of the field variable is missing
      outputFilePartNo++;

      outputFileParts[outputFilePartNo] << std::endl << std::string(4, '\t') << "</DataArray>";
    }

    outputFileParts[outputFilePartNo] << std::string(3, '\t') << "</PointData>" << std::endl
      << std::string(3, '\t') << "<CellData>" << std::endl
      << std::string(3, '\t') << "</CellData>" << std::endl
      << std::string(3, '\t') << "<Points>" << std::endl
      << std::string(4, '\t') << "<DataArray NumberOfComponents=\"3\" format=\"" << (binaryOutput_? "binary" : "ascii")
      << "\" >" << std::endl << std::string(5, '\t');

    // at this point the data of points (geometry field) is missing
    outputFilePartNo++;


    outputFileParts[outputFilePartNo] << std::string(4, '\t') << "</DataArray>" << std::endl
      << std::string(3, '\t') << "</Points>" << std::endl
      << std::string(3, '\t') << "<Verts></Verts>" << std::endl
      << std::string(3, '\t') << "<Lines>";

    if (vtkPieceIter->nVTKLines != 0)
    {
      // output the structural information of the lines
      // add connectivity array
      std::vector<int> values(2*vtkPieceIter->nVTKLines);
      for (int i = 0; i < vtkPieceIter->nVTKLines; i++)
      {
        values[2*i + 0] = i;
        values[2*i + 0] = i+1;
      }

      outputFileParts[outputFilePartNo] << std::endl << std::string(4, '\t') << "<DataArray Name=\"connectivity\" type=\"Int32\" ;

      if (binaryOutput_)
      {
        stringData = Paraview::encodeBase64(values);
        outputFileParts[outputFilePartNo] << "format=\"binary\" >" << std::endl;
      }
      else
      {
        stringData = Paraview::convertToAscii(values, fixedFormat);
        outputFileParts[outputFilePartNo] << "format=\"ascii\" >" << std::endl;
      }

      outputFileParts[outputFilePartNo] << std::string(5, '\t') << stringData << std::endl
        << std::string(4, '\t') << "</DataArray>" << std::endl;


      // add offsets array
      values.resize(vtkPieceIter->nVTKLines);
      for (int i = 0; i < vtkPieceIter->nVTKLines; i++)
      {
        values[i] = 2*i+1;
      }

      // output the structural information of the lines
      outputFileParts[outputFilePartNo] << std::endl << std::string(4, '\t') << "<DataArray Name=\"offsets\" type=\"Int32\" ;

      if (binaryOutput_)
      {
        stringData = Paraview::encodeBase64(values);
        outputFileParts[outputFilePartNo] << "format=\"binary\" >" << std::endl;
      }
      else
      {
        stringData = Paraview::convertToAscii(values, fixedFormat);
        outputFileParts[outputFilePartNo] << "format=\"ascii\" >" << std::endl;
      }

      outputFileParts[outputFilePartNo] << std::string(5, '\t') << stringData << std::endl
        << std::string(4, '\t') << "</DataArray>" << std::endl;
    }

    outputFileParts[outputFilePartNo] << "</Lines>" << std::endl
      << std::string(3, '\t') << "<Strips></Strips>" << std::endl
      << std::string(3, '\t') << "<Polys>";

    if (vtkPieceIter->nVTKPolys != 0)
    {
      // output the structural information of the polys
      // add connectivity array
      std::vector<int> values(2*vtkPieceIter->nVTKLines);
      for (int i = 0; i < vtkPieceIter->nVTKLines; i++)
      {
        values[2*i + 0] = i;
        values[2*i + 0] = i+1;
      }

      outputFileParts[outputFilePartNo] << std::endl << std::string(4, '\t') << "<DataArray Name=\"connectivity\" type=\"Int32\" ;

      if (binaryOutput_)
      {
        stringData = Paraview::encodeBase64(values);
        outputFileParts[outputFilePartNo] << "format=\"binary\" >" << std::endl;
      }
      else
      {
        stringData = Paraview::convertToAscii(values, fixedFormat);
        outputFileParts[outputFilePartNo] << "format=\"ascii\" >" << std::endl;
      }

      outputFileParts[outputFilePartNo] << std::string(5, '\t') << stringData << std::endl
        << std::string(4, '\t') << "</DataArray>" << std::endl;


      // add offsets array
      values.resize(vtkPieceIter->nVTKLines);
      for (int i = 0; i < vtkPieceIter->nVTKLines; i++)
      {
        values[i] = 2*i+1;
      }

      // output the structural information of the lines
      outputFileParts[outputFilePartNo] << std::endl << std::string(4, '\t') << "<DataArray Name=\"offsets\" type=\"Int32\" ;

      if (binaryOutput_)
      {
        stringData = Paraview::encodeBase64(values);
        outputFileParts[outputFilePartNo] << "format=\"binary\" >" << std::endl;
      }
      else
      {
        stringData = Paraview::convertToAscii(values, fixedFormat);
        outputFileParts[outputFilePartNo] << "format=\"ascii\" >" << std::endl;
      }

      outputFileParts[outputFilePartNo] << std::string(5, '\t') << stringData << std::endl
        << std::string(4, '\t') << "</DataArray>" << std::endl;
    }

    outputFileParts[outputFilePartNo] << "</Polys>" << std::endl
      << std::string(2, '\t') << "</Piece>" << std::endl;
  }

  outputFileParts[outputFilePartNo] << std::string(2, '\t') << "</Piece>" << std::endl
    << std::string(1, '\t') << "</PolyData>" << std::endl
    << "</VTKFile>" << std::endl;

  assert(outputFilePartNo+1 == nOutputFileParts);

  // loop over output file parts and collect the missing data for the own rank
}

};
