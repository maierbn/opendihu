#include "output_writer/paraview/paraview.h"

#include <fstream>
#include <iostream>
#include <thread>
#include <chrono>
#include <cstdio>  // remove

#include "easylogging++.h"
#include "base64.h"

#include "output_writer/paraview/loop_output.h"
#include "output_writer/paraview/loop_collect_mesh_properties.h"
#include "output_writer/paraview/loop_get_nodal_values.h"
#include "output_writer/paraview/loop_get_geometry_field_nodal_values.h"
#include "output_writer/paraview/poly_data_properties_for_mesh.h"
#include "control/diagnostic_tool/performance_measurement.h"

namespace OutputWriter
{

template<typename FieldVariablesForOutputWriterType>
void Paraview::writeCombinedUnstructuredGridFile(const FieldVariablesForOutputWriterType &fieldVariables, std::set<std::string> &meshNames,
                                                 bool output3DMeshes)
{
  // output a *.vtu file which contains 3D (if output3DMeshes==true) or 2D meshes (if output3DMeshes==false), if there are any

  std::map<std::string, PolyDataPropertiesForMesh> &meshPropertiesUnstructuredGridFile_ = (output3DMeshes? meshPropertiesUnstructuredGridFile3D_ : meshPropertiesUnstructuredGridFile2D_);

  bool meshPropertiesInitialized = !meshPropertiesUnstructuredGridFile_.empty();
  std::vector<std::string> meshNamesVector;

  // a counter that increments over every call to writeCombinedValuesVector per timestep and is used to access cached information
  int callIdentifier = 0;

  if (!meshPropertiesInitialized)
  {
    LOG(DEBUG) << "initialize meshPropertiesUnstructuredGridFile_";
    Control::PerformanceMeasurement::start("durationParaview3DInit");

    // collect the size data that is needed to compute offsets for parallel file output
    ParaviewLoopOverTuple::loopCollectMeshProperties<FieldVariablesForOutputWriterType>(fieldVariables, meshPropertiesUnstructuredGridFile_, meshNamesVector);

    Control::PerformanceMeasurement::stop("durationParaview3DInit");
  }
  else LOG(DEBUG) << "meshPropertiesUnstructuredGridFile_ already initialized";

  VLOG(1) << "writeCombinedUnstructuredGridFile on rankSubset_: " << *this->rankSubset_;
  assert(this->rankSubset_);

  VLOG(1) << "meshPropertiesUnstructuredGridFile_: " << meshPropertiesUnstructuredGridFile_;

  int targetDimensionality = 2;
  if (output3DMeshes)
    targetDimensionality = 3;

  VLOG(1) << "targetDimensionality: " << targetDimensionality;


  if (!meshPropertiesInitialized)
  {
    Control::PerformanceMeasurement::start("durationParaview2D3DInit");

    // parse the collected properties of the meshes that will be output to the file
    for (std::string meshName : meshNamesVector)
    {
      PolyDataPropertiesForMesh &polyDataPropertiesForMesh = meshPropertiesPolyDataFile_[meshName];

      // do not combine meshes other than 1D meshes
      if (polyDataPropertiesForMesh.dimensionality != targetDimensionality)
        continue;

      // check if this mesh should be combined with other meshes
      bool combineMesh = true;

      // check if mesh can be merged into previous meshes
      if (!vtkPiece3D_.properties.pointDataArrays.empty())   // if properties are already assigned by an earlier mesh
      {
        if (vtkPiece3D_.properties.pointDataArrays.size() != polyDataPropertiesForMesh.pointDataArrays.size())
        {
          LOG(DEBUG) << "Mesh " << meshName << " cannot be combined with " << vtkPiece3D_.meshNamesCombinedMeshes << ". Number of field variables mismatches for "
            << meshName << " (is " << polyDataPropertiesForMesh.pointDataArrays.size() << " instead of " << vtkPiece3D_.properties.pointDataArrays.size() << ")";
          combineMesh = false;
        }
        else
        {
          for (int j = 0; j < polyDataPropertiesForMesh.pointDataArrays.size(); j++)
          {
            if (vtkPiece3D_.properties.pointDataArrays[j].name != polyDataPropertiesForMesh.pointDataArrays[j].name)  // if the name of the jth field variable is different
            {
              LOG(DEBUG) << "Mesh " << meshName << " cannot be combined with " << vtkPiece3D_.meshNamesCombinedMeshes << ". Field variable names mismatch for "
                << meshName << " (there is \"" << vtkPiece3D_.properties.pointDataArrays[j].name << "\" instead of \"" << polyDataPropertiesForMesh.pointDataArrays[j].name << "\")";
              combineMesh = false;
            }
          }
        }

        if (combineMesh)
        {
          VLOG(1) << "Combine mesh " << meshName << " with " << vtkPiece3D_.meshNamesCombinedMeshes
            << ", add " << polyDataPropertiesForMesh.nPointsLocal << " points, " << polyDataPropertiesForMesh.nCellsLocal << " elements to "
            << vtkPiece3D_.properties.nPointsLocal << " points, " << vtkPiece3D_.properties.nCellsLocal << " elements";

          vtkPiece3D_.properties.nPointsLocal += polyDataPropertiesForMesh.nPointsLocal;
          vtkPiece3D_.properties.nCellsLocal += polyDataPropertiesForMesh.nCellsLocal;

          vtkPiece3D_.properties.nPointsGlobal += polyDataPropertiesForMesh.nPointsGlobal;
          vtkPiece3D_.properties.nCellsGlobal += polyDataPropertiesForMesh.nCellsGlobal;
          vtkPiece3D_.setVTKValues();
        }
      }
      else
      {
        VLOG(1) << "this is the first " << targetDimensionality << "D mesh";

        // properties are not yet assigned
        vtkPiece3D_.properties = polyDataPropertiesForMesh; // store properties
        vtkPiece3D_.setVTKValues();
      }

      VLOG(1) << "combineMesh: " << combineMesh;
      if (combineMesh)
      {
        // if the mesh is not yet present in meshNamesCombinedMeshes, add it
        if (vtkPiece3D_.meshNamesCombinedMeshes.find(meshName) == vtkPiece3D_.meshNamesCombinedMeshes.end())
        {
          vtkPiece3D_.meshNamesCombinedMeshes.insert(meshName);
          vtkPiece3D_.meshNamesCombinedMeshesVector.push_back(meshName);
        }
      }
    }

    LOG(DEBUG) << "vtkPiece3D_: meshNamesCombinedMeshes: " << vtkPiece3D_.meshNamesCombinedMeshes << ", properties: " << vtkPiece3D_.properties
      << ", firstScalarName: " << vtkPiece3D_.firstScalarName << ", firstVectorName: " << vtkPiece3D_.firstVectorName;

    Control::PerformanceMeasurement::stop("durationParaview2D3DInit");
  }

  // write combined mesh
  if (!vtkPiece3D_.meshNamesCombinedMeshes.empty())
  {
    PolyDataPropertiesForMesh &polyDataPropertiesForMesh = vtkPiece3D_.properties;

    VLOG(1) << "polyDataPropertiesForMesh combined: " << polyDataPropertiesForMesh;

    // if the meshes have the right dimensionality (2D or 3D)
    if (vtkPiece3D_.properties.dimensionality == targetDimensionality)
    {
      meshNames = vtkPiece3D_.meshNamesCombinedMeshes;

      std::stringstream filename;
      filename << this->filenameBaseWithNo_ << ".vtu";

      // write actual file
      writeCombinedUnstructuredGridFile<FieldVariablesForOutputWriterType>(
        fieldVariables, vtkPiece3D_.properties, meshPropertiesUnstructuredGridFile_, vtkPiece3D_.meshNamesCombinedMeshesVector, meshPropertiesInitialized, callIdentifier, filename.str());
      callIdentifier++;
    }
    else
    {
      VLOG(1) << "skip meshes " << vtkPiece3D_.meshNamesCombinedMeshes << " because " << vtkPiece3D_.properties.dimensionality << " != " << targetDimensionality;
    }
  }

  // write other single meshes

  // loop over 3D or 2D meshes
  for (std::map<std::string, PolyDataPropertiesForMesh>::iterator meshPropertiesIter = meshPropertiesUnstructuredGridFile_.begin(); meshPropertiesIter != meshPropertiesUnstructuredGridFile_.end(); meshPropertiesIter++)
  {
    PolyDataPropertiesForMesh &polyDataPropertiesForMesh = meshPropertiesIter->second;

    // if the current mesh does not have the correct dimensionality, skip this mesh
    if (polyDataPropertiesForMesh.dimensionality != targetDimensionality)
      continue;

    // if the current mesh was already output within the combined file, skip this mesh
    if (meshNames.find(meshPropertiesIter->first) != meshNames.end())
      continue;

    std::stringstream filename;
    if (meshPropertiesUnstructuredGridFile_.size() > 1)
    {
      // split the current filename into base part and number part, e.g. "output_0000010" -> "output" + "_0000010"
      std::string filenamePartBase;
      std::string filenamePartNumber;
      
      if (this->filenameBaseWithNo_.find("_") != std::string::npos)
      {
        filenamePartBase = this->filenameBaseWithNo_.substr(0, this->filenameBaseWithNo_.rfind("_"));
        filenamePartNumber = this->filenameBaseWithNo_.substr(this->filenameBaseWithNo_.rfind("_"));
      }
      else 
      {
        filenamePartBase = this->filenameBaseWithNo_;
      }

      // add mesh name in filename
      filename << filenamePartBase << "_" << meshPropertiesIter->first << filenamePartNumber << ".vtu";
    }
    else
    {
      std::stringstream filename;
      filename << this->filenameBaseWithNo_ << ".vtu";
    }

    LOG(DEBUG) << "meshNames: " << meshNames << ", next mesh to write: " << meshPropertiesIter->first << ", filename: " << filename.str();

    // write actual file
    std::vector<std::string> currentMesh;
    currentMesh.push_back(meshPropertiesIter->first);
    meshNames.insert(meshPropertiesIter->first);

    writeCombinedUnstructuredGridFile<FieldVariablesForOutputWriterType>(
      fieldVariables, polyDataPropertiesForMesh, meshPropertiesUnstructuredGridFile_, currentMesh, meshPropertiesInitialized, callIdentifier, filename.str());
    callIdentifier++;

    // for the next meshes, reinitialize, the initialized information can only be shared if there is exactly 1 mesh to write in this method
    meshPropertiesInitialized = false;
  }
}


template<typename FieldVariablesForOutputWriterType>
void Paraview::writeCombinedUnstructuredGridFile(const FieldVariablesForOutputWriterType &fieldVariables, PolyDataPropertiesForMesh &polyDataPropertiesForMesh,
                                                 const std::map<std::string, PolyDataPropertiesForMesh> &meshPropertiesUnstructuredGridFile,
                                                 std::vector<std::string> meshNames,
                                                 bool meshPropertiesInitialized, int &callIdentifier, std::string filename)
{
  int targetDimensionality = polyDataPropertiesForMesh.dimensionality;
  bool output3DMeshes = targetDimensionality == 3;

  std::set<std::string> meshNamesSet(meshNames.begin(), meshNames.end());

  VLOG(1) << "writeCombinedUnstructuredGridFile, filename=" << filename << ", meshNames: " << meshNames << ", meshPropertiesInitialized=" << meshPropertiesInitialized << ", targetDimensionality: " << targetDimensionality;

  //! find out name of first scalar and vector field variables
  std::string firstScalarName;
  std::string firstVectorName;

  for (std::vector<PolyDataPropertiesForMesh::DataArrayName>::iterator iter = polyDataPropertiesForMesh.pointDataArrays.begin(); iter != polyDataPropertiesForMesh.pointDataArrays.end(); iter++)
  {
    if (iter->nComponents == 3 && firstVectorName.empty())
    {
      firstVectorName = iter->name;
    }
    else if (firstScalarName.empty())
    {
      firstScalarName = iter->name;
    }

    if (!firstScalarName.empty() && firstVectorName.empty())
      break;
  }

  // determine filename, broadcast from rank 0
  int filenameLength = filename.length();

  // broadcast length of filename
  MPIUtility::handleReturnValue(MPI_Bcast(&filenameLength, 1, MPI_INT, 0, this->rankSubset_->mpiCommunicator()), "MPI_Bcast (5)");

  std::vector<char> receiveBuffer(filenameLength+1, char(0));
  strcpy(receiveBuffer.data(), filename.c_str());
  MPIUtility::handleReturnValue(MPI_Bcast(receiveBuffer.data(), filenameLength, MPI_CHAR, 0, this->rankSubset_->mpiCommunicator()), "MPI_Bcast (6)");

  std::string filenameStr(receiveBuffer.begin(), receiveBuffer.end());

  // remove file if it exists, synchronization afterwards by MPI calls, that is why the remove call is already here
  assert(this->rankSubset_);
  int ownRankNo = this->rankSubset_->ownRankNo();
  if (ownRankNo == 0)
  {
    // open file to ensure that directory exists and file is writable
    std::ofstream file;
    Generic::openFile(file, filenameStr);

    // close and delete file
    file.close();
    std::remove(filenameStr.c_str());
  }

  // exchange information about offset in terms of nCells and nPoints
  if (!meshPropertiesInitialized)
  {
    VLOG(1) << "set nCellsPreviousRanks3D_, nPointsPreviousRanks3D_, nPointsGlobal3D_";

    Control::PerformanceMeasurement::start("durationParaview3DInit");
    Control::PerformanceMeasurement::start("durationParaview3DReduction");
    nPointsGlobal3D_[callIdentifier] = 0;
    MPIUtility::handleReturnValue(MPI_Reduce(&polyDataPropertiesForMesh.nPointsLocal, &nPointsGlobal3D_[callIdentifier], 1, MPI_INT, MPI_SUM, 0, this->rankSubset_->mpiCommunicator()), "MPI_Reduce");

    for (std::string meshName : meshNames)
    {
      nCellsPreviousRanks3D_[meshName] = 0;
      nPointsPreviousRanks3D_[meshName] = 0;

      MPIUtility::handleReturnValue(MPI_Exscan(&meshPropertiesUnstructuredGridFile.at(meshName).nCellsLocal, &nCellsPreviousRanks3D_[meshName], 1, MPI_INT, MPI_SUM, this->rankSubset_->mpiCommunicator()), "MPI_Exscan");
      MPIUtility::handleReturnValue(MPI_Exscan(&meshPropertiesUnstructuredGridFile.at(meshName).nPointsLocal, &nPointsPreviousRanks3D_[meshName], 1, MPI_INT, MPI_SUM, this->rankSubset_->mpiCommunicator()), "MPI_Exscan");

      VLOG(1) << "meshName \"" << meshName << "\"";
      VLOG(1) << "nCellsLocal local: " << meshPropertiesUnstructuredGridFile.at(meshName).nCellsLocal << ", prefix sum: " << nCellsPreviousRanks3D_[meshName];
      VLOG(1) << "nPointsLocal local: " << meshPropertiesUnstructuredGridFile.at(meshName).nPointsLocal << ", prefix sum: " << nPointsPreviousRanks3D_[meshName];
    }

    Control::PerformanceMeasurement::stop("durationParaview3DReduction");
    Control::PerformanceMeasurement::stop("durationParaview3DInit");

    VLOG(1) << "nPointsGlobal: " << nPointsGlobal3D_[callIdentifier];
  }
  else
  {
    VLOG(1) << "nCellsPreviousRanks3D_, nPointsPreviousRanks3D_, nPointsGlobal3D_ already set to "
      << nCellsPreviousRanks3D_ << ", " << nPointsPreviousRanks3D_ << ", " << nPointsGlobal3D_;
  }

#ifndef NDEBUG
  std::vector<node_no_t> &nNodesLocalWithGhosts = polyDataPropertiesForMesh.nNodesLocalWithGhosts;
  assert(nNodesLocalWithGhosts.size() == targetDimensionality);
#endif

  int nNodesPerCell = 4;
  if (output3DMeshes)
    nNodesPerCell = 8;

  // get local data values
  // setup connectivity array, which gives the node numbers for every element/cell
  int nConnectivityValues = 0;

  // loop over the meshes that will be combined into the current file
  // meshPropertiesUnstructuredGridFile[...] contains information for all meshes, polyDataPropertiesForMesh contains already combined information of the set of meshes
  for (std::string meshName : meshNames)
  {
    nConnectivityValues += meshPropertiesUnstructuredGridFile.at(meshName).nCellsLocal * nNodesPerCell;
    VLOG(1) << "  mesh \"" << meshName << "\" has " << meshPropertiesUnstructuredGridFile.at(meshName).nCellsLocal << " local cells, " << nNodesPerCell
      << ", " << meshPropertiesUnstructuredGridFile.at(meshName).nCellsLocal * nNodesPerCell << " connectivity values for this mesh";
  }
  VLOG(1) << "nConnectivityValues: " << nConnectivityValues;

  std::vector<int> connectivityValues(nConnectivityValues);

  VLOG(1) << "n connectivity values from unstructured: " << polyDataPropertiesForMesh.unstructuredMeshConnectivityValues.size();
  VLOG(1) << "nCellsLocal: " << polyDataPropertiesForMesh.nCellsLocal;

  // if connectivity values are already explicitly given, this is the case if we have an unstructured mesh to output
  if (!polyDataPropertiesForMesh.unstructuredMeshConnectivityValues.empty())
  {
    if (meshNames.size() > 1)
      LOG(FATAL) << "Cannot combine file for unstructured meshes.";
    assert(polyDataPropertiesForMesh.unstructuredMeshConnectivityValues.size() == connectivityValues.size());

    VLOG(1) << "connectivityValues is initialized to " << connectivityValues.size() << ", values: " << connectivityValues;
    VLOG(1) << "now copy " << polyDataPropertiesForMesh.unstructuredMeshConnectivityValues.size() << ", values: " << polyDataPropertiesForMesh.unstructuredMeshConnectivityValues;
    std::copy(polyDataPropertiesForMesh.unstructuredMeshConnectivityValues.begin(), polyDataPropertiesForMesh.unstructuredMeshConnectivityValues.end(), connectivityValues.begin());
  }
  else
  {
    int connectivityValuesIndexOffset = 0;
    int connectivityValuesOffset = 0;

    for (std::string meshName : meshNames)
    {
      connectivityValuesOffset += nPointsPreviousRanks3D_[meshName];
    }

    VLOG(1) << "connectivityValuesOffset start: " << connectivityValuesOffset;

    // loop over meshes
    for (std::string meshName : meshNames)
    {
      int nPointsPreviousRanks3D = connectivityValuesOffset;

      const std::vector<node_no_t> &nNodesLocalWithGhosts = meshPropertiesUnstructuredGridFile.at(meshName).nNodesLocalWithGhosts;

      VLOG(1) << "connectivity for mesh \"" << meshName << "\",  nPointsPreviousRanks3D: " << nPointsPreviousRanks3D << ", nNodesLocalWithGhosts: " << nNodesLocalWithGhosts;
      VLOG(1) << "connectivityValuesIndexOffset: " << connectivityValuesIndexOffset;
      VLOG(1) << "nNodesLocalWithGhosts: " << nNodesLocalWithGhosts;

      // for structured meshes create connectivity values now
      if (output3DMeshes)
      {
        // fill connectivityValues for 3D meshes
        element_no_t elementIndex = 0;
        for (int indexZ = 0; indexZ < nNodesLocalWithGhosts[2]-1; indexZ++)
        {
          for (int indexY = 0; indexY < nNodesLocalWithGhosts[1]-1; indexY++)
          {
            for (int indexX = 0; indexX < nNodesLocalWithGhosts[0]-1; indexX++, elementIndex++)
            {
              if (connectivityValuesIndexOffset + elementIndex*8 + 7 >= connectivityValues.size())
              {
                LOG(FATAL) << connectivityValuesIndexOffset + elementIndex*8 + 7 << ">= " << connectivityValues.size()
                << ", connectivityValues are not large enough: " << connectivityValues.size() << ", but "
                  << connectivityValuesIndexOffset << "+" << nNodesLocalWithGhosts[0]-1 << "x" << nNodesLocalWithGhosts[1]-1 << "x" << nNodesLocalWithGhosts[2]-1 << " = "
                  << connectivityValuesIndexOffset + (nNodesLocalWithGhosts[0]-1)*(nNodesLocalWithGhosts[1]-1)*(nNodesLocalWithGhosts[2]-1) << " elements";
              }

              connectivityValues[connectivityValuesIndexOffset + elementIndex*8 + 0]
                = nPointsPreviousRanks3D
                + indexZ*nNodesLocalWithGhosts[0]*nNodesLocalWithGhosts[1]
                + indexY*nNodesLocalWithGhosts[0] + indexX;
              connectivityValues[connectivityValuesIndexOffset + elementIndex*8 + 1]
                = nPointsPreviousRanks3D
                + indexZ*nNodesLocalWithGhosts[0]*nNodesLocalWithGhosts[1]
                + indexY*nNodesLocalWithGhosts[0] + indexX + 1;
              connectivityValues[connectivityValuesIndexOffset + elementIndex*8 + 2]
                = nPointsPreviousRanks3D
                + indexZ*nNodesLocalWithGhosts[0]*nNodesLocalWithGhosts[1]
                + (indexY+1)*nNodesLocalWithGhosts[0] + indexX + 1;
              connectivityValues[connectivityValuesIndexOffset + elementIndex*8 + 3]
                = nPointsPreviousRanks3D
                + indexZ*nNodesLocalWithGhosts[0]*nNodesLocalWithGhosts[1]
                + (indexY+1)*nNodesLocalWithGhosts[0] + indexX;
              connectivityValues[connectivityValuesIndexOffset + elementIndex*8 + 4]
                = nPointsPreviousRanks3D
                + (indexZ+1)*nNodesLocalWithGhosts[0]*nNodesLocalWithGhosts[1]
                + indexY*nNodesLocalWithGhosts[0] + indexX;
              connectivityValues[connectivityValuesIndexOffset + elementIndex*8 + 5]
                = nPointsPreviousRanks3D
                + (indexZ+1)*nNodesLocalWithGhosts[0]*nNodesLocalWithGhosts[1]
                + indexY*nNodesLocalWithGhosts[0] + indexX + 1;
              connectivityValues[connectivityValuesIndexOffset + elementIndex*8 + 6]
                = nPointsPreviousRanks3D
                + (indexZ+1)*nNodesLocalWithGhosts[0]*nNodesLocalWithGhosts[1]
                + (indexY+1)*nNodesLocalWithGhosts[0] + indexX + 1;
              connectivityValues[connectivityValuesIndexOffset + elementIndex*8 + 7]
                = nPointsPreviousRanks3D
                + (indexZ+1)*nNodesLocalWithGhosts[0]*nNodesLocalWithGhosts[1]
                + (indexY+1)*nNodesLocalWithGhosts[0] + indexX;
            }
          }
        }
      }
      else
      {
        // fill connectivityValues for 2D meshes
        element_no_t elementIndex = 0;
        for (int indexY = 0; indexY < nNodesLocalWithGhosts[1]-1; indexY++)
        {
          for (int indexX = 0; indexX < nNodesLocalWithGhosts[0]-1; indexX++, elementIndex++)
          {
            if (connectivityValuesIndexOffset + elementIndex*4 + 3 >= connectivityValues.size())
            {
              LOG(FATAL) << connectivityValuesIndexOffset + elementIndex*4 + 3 << ">= " << connectivityValues.size()
                << ", connectivityValues are not large enough: " << connectivityValues.size() << ", but "
                << connectivityValuesIndexOffset << " + " << nNodesLocalWithGhosts[0]-1 << "x" << nNodesLocalWithGhosts[1]-1 << " = "
                << connectivityValuesIndexOffset + (nNodesLocalWithGhosts[0]-1)*(nNodesLocalWithGhosts[1]-1) << " elements";
            }

            connectivityValues[connectivityValuesIndexOffset + elementIndex*4 + 0]
              = nPointsPreviousRanks3D
              + indexY*nNodesLocalWithGhosts[0] + indexX;
            connectivityValues[connectivityValuesIndexOffset + elementIndex*4 + 1]
              = nPointsPreviousRanks3D
              + indexY*nNodesLocalWithGhosts[0] + indexX + 1;
            connectivityValues[connectivityValuesIndexOffset + elementIndex*4 + 2]
              = nPointsPreviousRanks3D
              + (indexY+1)*nNodesLocalWithGhosts[0] + indexX + 1;
            connectivityValues[connectivityValuesIndexOffset + elementIndex*4 + 3]
              = nPointsPreviousRanks3D
              + (indexY+1)*nNodesLocalWithGhosts[0] + indexX;
          }
        }
      }

      connectivityValuesIndexOffset += meshPropertiesUnstructuredGridFile.at(meshName).nCellsLocal * nNodesPerCell;
      connectivityValuesOffset += meshPropertiesUnstructuredGridFile.at(meshName).nPointsLocal;
    }
  }

  VLOG(1) << "nPointsPreviousRanks3D_: " << nPointsPreviousRanks3D_ << ", nCellsPreviousRanks3D_: " << nCellsPreviousRanks3D_;
  //VLOG(1) << "nNodesLocalWithGhosts: " << nNodesLocalWithGhosts[0]-1 << "x" << nNodesLocalWithGhosts[1]-1 << "x" << nNodesLocalWithGhosts[2]-1;
  VLOG(1) << "connectivity: " << connectivityValues;

  // setup offset array
  std::vector<int> offsetValues(polyDataPropertiesForMesh.nCellsLocal);

  int offsetValuesIndexOffset = 0;
  int offsetValuesOffset = 0;

  for (std::string meshName : meshNames)
  {
    offsetValuesOffset += nCellsPreviousRanks3D_[meshName] * nNodesPerCell;
  }

  // loop over meshes
  for (std::string meshName : meshNames)
  {
    for (int i = 0; i < meshPropertiesUnstructuredGridFile.at(meshName).nCellsLocal; i++)
    {
      assert(offsetValuesIndexOffset + i < offsetValues.size());
      offsetValues[offsetValuesIndexOffset + i] = offsetValuesOffset + (i + 1) * nNodesPerCell;    // specifies the end, i.e. one after the last, of the last of nodes for each element
    }

    offsetValuesIndexOffset += meshPropertiesUnstructuredGridFile.at(meshName).nCellsLocal;
    offsetValuesOffset += meshPropertiesUnstructuredGridFile.at(meshName).nCellsLocal * nNodesPerCell;
  }

  VLOG(1) << "offsetValues: " << offsetValues;

  // collect all data for the field variables, organized by field variable names
  std::map<std::string, std::vector<double>> fieldVariableValues;
  ParaviewLoopOverTuple::loopGetNodalValues<FieldVariablesForOutputWriterType>(fieldVariables, meshNamesSet, fieldVariableValues);

  if (!meshPropertiesInitialized)
  {
    // if next assertion fails, output why for debugging
    if (fieldVariableValues.size() != polyDataPropertiesForMesh.pointDataArrays.size())
    {
      LOG(DEBUG) << "n field variable values: " << fieldVariableValues.size() << ", n point data arrays: "
        << polyDataPropertiesForMesh.pointDataArrays.size();
      LOG(DEBUG) << "mesh names: " << meshNames;
      std::stringstream pointDataArraysNames;
      for (int i = 0; i < polyDataPropertiesForMesh.pointDataArrays.size(); i++)
      {
        pointDataArraysNames << polyDataPropertiesForMesh.pointDataArrays[i].name << " ";
      }
      LOG(DEBUG) << "pointDataArraysNames: " <<  pointDataArraysNames.str();
      LOG(DEBUG) << "FieldVariablesForOutputWriterType: " << StringUtility::demangle(typeid(FieldVariablesForOutputWriterType).name());
    }

    assert(fieldVariableValues.size() == polyDataPropertiesForMesh.pointDataArrays.size());
  }

  // output 3D or 2D mesh
  if (!meshPropertiesInitialized)
  {
    // add field variable "partitioning" with 1 component
    PolyDataPropertiesForMesh::DataArrayName dataArrayName;
    dataArrayName.name = "partitioning";
    dataArrayName.nComponents = 1;
    dataArrayName.componentNames = std::vector<std::string>(1, "rankNo");

    polyDataPropertiesForMesh.pointDataArrays.push_back(dataArrayName);
  }

  // set data for partitioning field variable
  assert (!fieldVariableValues.empty());
  fieldVariableValues["partitioning"].resize(polyDataPropertiesForMesh.nPointsLocal, (double)this->rankSubset_->ownRankNo());

  // for 2D field variables, add 0 every 2nd entry to generate 3D values, because paraview cannot handle 2D vectors properly
  if (targetDimensionality == 2)
  {
    std::vector<double> buffer;

    // loop over all field variables
    for (std::vector<PolyDataPropertiesForMesh::DataArrayName>::iterator pointDataArrayIter = polyDataPropertiesForMesh.pointDataArrays.begin(); pointDataArrayIter != polyDataPropertiesForMesh.pointDataArrays.end(); pointDataArrayIter++)
    {
      // if it is a 2D vector field
      if (pointDataArrayIter->nComponents == 2)
      {
        std::vector<double> &values = fieldVariableValues[pointDataArrayIter->name];

        // copy all values to a buffer
        buffer.assign(values.begin(), values.end());
        int nValues = buffer.size();

        // resize the values vector to 2/3 the size
        values.resize(nValues/2*3);

        // loop over the new values vector and set the entries
        for (int i = 0; i < nValues/2; i++)
        {
          values[3*i + 0] = buffer[2*i + 0];
          values[3*i + 1] = buffer[2*i + 1];
          values[3*i + 2] = 0.0;
        }
      }
    }
  }

  // collect all data for the geometry field variable
  std::vector<double> geometryFieldValues;
  ParaviewLoopOverTuple::loopGetGeometryFieldNodalValues<FieldVariablesForOutputWriterType>(fieldVariables, meshNamesSet, geometryFieldValues);

  VLOG(1) << "meshNames: " << meshNames << ", rank " << this->rankSubset_->ownRankNo() << ", n geometryFieldValues: " << geometryFieldValues.size();
  if (geometryFieldValues.size() == 0)
  {
    LOG(FATAL) << "There is no geometry field. You have to provide a geomteryField in the field variables returned by getFieldVariablesForOutputWriter!";
  }

  int nOutputFileParts = 5 + polyDataPropertiesForMesh.pointDataArrays.size();

  // transform current time to string
  std::vector<double> time(1, this->currentTime_);
  std::string stringTime;
  if (binaryOutput_)
  {
    stringTime = Paraview::encodeBase64Float(time.begin(), time.end());
  }
  else
  {
    stringTime = Paraview::convertToAscii(time, fixedFormat_);
  }

  // create the basic structure of the output file
  std::vector<std::stringstream> outputFileParts(nOutputFileParts);
  int outputFilePartNo = 0;
  outputFileParts[outputFilePartNo] << "<?xml version=\"1.0\"?>" << std::endl
    << "<!-- " << DihuContext::versionText() << " " << DihuContext::metaText()
    << ", currentTime: " << this->currentTime_ << ", timeStepNo: " << this->timeStepNo_ << " -->" << std::endl
    << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">" << std::endl    // intel cpus are LittleEndian
    << std::string(1, '\t') << "<UnstructuredGrid>" << std::endl
    << std::string(2, '\t') << "<FieldData>" << std::endl
    << std::string(3, '\t') << "<DataArray type=\"Float32\" Name=\"Time\" NumberOfTuples=\"1\" format=\"" << (binaryOutput_? "binary" : "ascii")
    << "\" >" << std::endl
    << std::string(4, '\t') << stringTime << std::endl
    << std::string(3, '\t') << "</DataArray>" << std::endl
    << std::string(2, '\t') << "</FieldData>" << std::endl;

  outputFileParts[outputFilePartNo] << std::string(2, '\t') << "<Piece NumberOfPoints=\"" << nPointsGlobal3D_[callIdentifier]
    << "\" NumberOfCells=\"" << polyDataPropertiesForMesh.nCellsGlobal << "\">" << std::endl
    << std::string(3, '\t') << "<PointData";

  if (firstScalarName != "")
  {
    outputFileParts[outputFilePartNo] << " Scalars=\"" << firstScalarName << "\"";
  }
  if (firstVectorName != "")
  {
    outputFileParts[outputFilePartNo] << " Vectors=\"" << firstVectorName << "\"";
  }
  outputFileParts[outputFilePartNo] << ">" << std::endl;

  // loop over field variables (PointData)
  for (std::vector<PolyDataPropertiesForMesh::DataArrayName>::iterator pointDataArrayIter = polyDataPropertiesForMesh.pointDataArrays.begin(); pointDataArrayIter != polyDataPropertiesForMesh.pointDataArrays.end(); pointDataArrayIter++)
  {
    // paraview cannot handle 2D vector fields without warnings, so set 3D vector fields
    int nComponentsParaview = pointDataArrayIter->nComponents;
    if (nComponentsParaview == 2)
      nComponentsParaview = 3;

    // set up string for component names
    std::stringstream componentNames;
    bool isComponentNamesSet = false;
    for (int componentNo = 0; componentNo < pointDataArrayIter->nComponents; componentNo++)
    {
      // get component name of the current component
      std::string componentName = pointDataArrayIter->componentNames[componentNo];

      componentNames << "ComponentName" << componentNo << "=\"" << componentName << "\" ";

      // check if it is equal to simply the componentNo, if no component name was explicitly specify they have names "0", "1", ...
      std::stringstream trivialComponentName;
      trivialComponentName << componentNo;

      if (componentName != trivialComponentName.str())
        isComponentNamesSet = true;
    }

    // if there were no real component names set, do not include them in the VTK file
    if (!isComponentNamesSet)
      componentNames.str("");

    // write normal data element
    outputFileParts[outputFilePartNo] << std::string(4, '\t') << "<DataArray "
        << "Name=\"" << pointDataArrayIter->name << "\" "
        << "type=\"" << (pointDataArrayIter->name == "partitioning"? "Int32" : "Float32") << "\" "
        << "NumberOfComponents=\"" << nComponentsParaview << "\" "
        << componentNames.str()
        << " format=\"" << (binaryOutput_? "binary" : "ascii")
        << "\" >" << std::endl << std::string(5, '\t');

    // at this point the data of the field variable is missing
    outputFilePartNo++;

    outputFileParts[outputFilePartNo] << std::endl << std::string(4, '\t') << "</DataArray>" << std::endl;
  }

  outputFileParts[outputFilePartNo] << std::string(3, '\t') << "</PointData>" << std::endl
    << std::string(3, '\t') << "<CellData>" << std::endl
    << std::string(3, '\t') << "</CellData>" << std::endl
    << std::string(3, '\t') << "<Points>" << std::endl
    << std::string(4, '\t') << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"" << (binaryOutput_? "binary" : "ascii")
    << "\" >" << std::endl << std::string(5, '\t');

  // at this point the data of points (geometry field) is missing
  outputFilePartNo++;

  outputFileParts[outputFilePartNo] << std::endl << std::string(4, '\t') << "</DataArray>" << std::endl
    << std::string(3, '\t') << "</Points>" << std::endl
    << std::string(3, '\t') << "<Cells>" << std::endl
    << std::string(4, '\t') << "<DataArray Name=\"connectivity\" type=\"Int32\" "
    << (binaryOutput_? "format=\"binary\"" : "format=\"ascii\"") << ">" << std::endl << std::string(5, '\t');

  // at this point the the structural information of the lines (connectivity) is missing
  outputFilePartNo++;

  outputFileParts[outputFilePartNo]
    << std::endl << std::string(4, '\t') << "</DataArray>" << std::endl
    << std::string(4, '\t') << "<DataArray Name=\"offsets\" type=\"Int32\" "
    << (binaryOutput_? "format=\"binary\"" : "format=\"ascii\"") << ">" << std::endl << std::string(5, '\t');

  // at this point the offset array will be written to the file
  outputFilePartNo++;

  outputFileParts[outputFilePartNo]
    << std::endl << std::string(4, '\t') << "</DataArray>" << std::endl
    << std::string(4, '\t') << "<DataArray Name=\"types\" type=\"UInt8\" "
    << (binaryOutput_? "format=\"binary\"" : "format=\"ascii\"") << ">" << std::endl << std::string(5, '\t');

  // at this point the types array will be written to the file
  outputFilePartNo++;

  outputFileParts[outputFilePartNo]
    << std::endl << std::string(4, '\t') << "</DataArray>" << std::endl
    << std::string(3, '\t') << "</Cells>" << std::endl
    << std::string(2, '\t') << "</Piece>" << std::endl
    << std::string(1, '\t') << "</UnstructuredGrid>" << std::endl
    << "</VTKFile>" << std::endl;

  assert(outputFilePartNo+1 == nOutputFileParts);

  // loop over output file parts and collect the missing data for the own rank

  VLOG(1) << "outputFileParts:";

  for (std::vector<std::stringstream>::iterator iter = outputFileParts.begin(); iter != outputFileParts.end(); iter++)
  {
    VLOG(1) << "  " << iter->str();
  }

  LOG(DEBUG) << "open MPI file \"" << filenameStr << "\" for rankSubset " << *this->rankSubset_;

  // open file
  MPI_File fileHandle;
  MPIUtility::handleReturnValue(MPI_File_open(this->rankSubset_->mpiCommunicator(), filenameStr.c_str(),
                                              //MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_UNIQUE_OPEN,
                                              MPI_MODE_WRONLY | MPI_MODE_CREATE,
                                              MPI_INFO_NULL, &fileHandle), "MPI_File_open");

  Control::PerformanceMeasurement::start("durationParaview3DWrite");

  // write beginning of file on rank 0
  outputFilePartNo = 0;

  writeAsciiDataShared(fileHandle, ownRankNo, outputFileParts[outputFilePartNo].str());
  outputFilePartNo++;

  VLOG(1) << "get current shared file position";

  // get current file position
  MPI_Offset currentFilePosition = 0;
  MPIUtility::handleReturnValue(MPI_File_get_position_shared(fileHandle, &currentFilePosition), "MPI_File_get_position_shared");
  LOG(DEBUG) << "current shared file position: " << currentFilePosition;

  // write field variables
  // loop over field variables
  for (std::vector<PolyDataPropertiesForMesh::DataArrayName>::iterator pointDataArrayIter = polyDataPropertiesForMesh.pointDataArrays.begin();
      pointDataArrayIter != polyDataPropertiesForMesh.pointDataArrays.end(); pointDataArrayIter++)
  {
    assert(fieldVariableValues.find(pointDataArrayIter->name) != fieldVariableValues.end());

    VLOG(1) << "write vector for field variable \"" << pointDataArrayIter->name << "\".";

    // write values
    bool writeFloatsAsInt = pointDataArrayIter->name == "partitioning";    // for partitioning, convert float values to integer values for output
    writeCombinedValuesVector(fileHandle, ownRankNo, fieldVariableValues[pointDataArrayIter->name], callIdentifier++, writeFloatsAsInt);

    // write next xml constructs
    writeAsciiDataShared(fileHandle, ownRankNo, outputFileParts[outputFilePartNo].str());
    outputFilePartNo++;
  }

  VLOG(1) << "write vector for geometry data";

  // write geometry field data
  writeCombinedValuesVector(fileHandle, ownRankNo, geometryFieldValues, callIdentifier++);

  // write next xml constructs
  writeAsciiDataShared(fileHandle, ownRankNo, outputFileParts[outputFilePartNo].str());
  outputFilePartNo++;

  // write connectivity values
  writeCombinedValuesVector(fileHandle, ownRankNo, connectivityValues, callIdentifier++);

  // write next xml constructs
  writeAsciiDataShared(fileHandle, ownRankNo, outputFileParts[outputFilePartNo].str());
  outputFilePartNo++;

  // write offset values
  writeCombinedValuesVector(fileHandle, ownRankNo, offsetValues, callIdentifier++);

  // write next xml constructs
  writeAsciiDataShared(fileHandle, ownRankNo, outputFileParts[outputFilePartNo].str());
  outputFilePartNo++;

  // write types values
  writeCombinedTypesVector(fileHandle, ownRankNo, polyDataPropertiesForMesh.nCellsGlobal, output3DMeshes, callIdentifier++);

  // write next xml constructs
  writeAsciiDataShared(fileHandle, ownRankNo, outputFileParts[outputFilePartNo].str());

  /*
    int array_of_sizes[1];
    array_of_sizes[0]=numProcs;
    int array_of_subsizes[1];
    array_of_subsizes[0]=1;
    int array_of_starts[1];
    array_of_starts[0]=myId;


    MPI_Datatype accessPattern;
    MPI_Type_create_subarray(1,array_of_sizes, array_of_subsizes, array_of_starts, MPI_ORDER_C, MPI_BYTE, &accessPattern);
    MPI_Type_commit(&accessPattern);

    MPI_File_set_view(fh, 0, MPI_BYTE, accessPattern, "native", MPI_INFO_NULL);
    MPI_File_write(fh, v, size, MPI_BYTE, MPI_STATUS_IGNORE);
  */

  Control::PerformanceMeasurement::stop("durationParaview3DWrite");
  MPIUtility::handleReturnValue(MPI_File_close(&fileHandle), "MPI_File_close");

  // register file at SeriesWriter to be included in the "*.vtk.series" JSON file
  if (ownRankNo == 0)
  {
    Paraview::seriesWriter().registerNewFile(filenameStr, this->currentTime_);
  }
}



} // namespace
