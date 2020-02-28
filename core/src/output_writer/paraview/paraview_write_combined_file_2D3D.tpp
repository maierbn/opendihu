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

  if (!meshPropertiesInitialized)
  {
    LOG(DEBUG) << "initialize meshPropertiesUnstructuredGridFile_";
    Control::PerformanceMeasurement::start("durationParaview3DInit");

    // collect the size data that is needed to compute offsets for parallel file output
    ParaviewLoopOverTuple::loopCollectMeshProperties<FieldVariablesForOutputWriterType>(fieldVariables, meshPropertiesUnstructuredGridFile_);

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
    for (std::map<std::string,PolyDataPropertiesForMesh>::iterator iter = meshPropertiesPolyDataFile_.begin(); iter != meshPropertiesPolyDataFile_.end(); iter++)
    {
      std::string meshName = iter->first;

      // do not combine meshes other than 1D meshes
      if (iter->second.dimensionality != targetDimensionality)
        continue;

      // check if this mesh should be combined with other meshes
      bool combineMesh = true;

      // check if mesh can be merged into previous meshes
      if (!vtkPiece_.properties.pointDataArrays.empty())   // if properties are already assigned by an earlier mesh
      {
        if (vtkPiece_.properties.pointDataArrays.size() != iter->second.pointDataArrays.size())
        {
          LOG(DEBUG) << "Mesh " << meshName << " cannot be combined with " << vtkPiece_.meshNamesCombinedMeshes << ". Number of field variables mismatches for "
            << meshName << " (is " << iter->second.pointDataArrays.size() << " instead of " << vtkPiece_.properties.pointDataArrays.size() << ")";
          combineMesh = false;
        }
        else
        {
          for (int j = 0; j < iter->second.pointDataArrays.size(); j++)
          {
            if (vtkPiece_.properties.pointDataArrays[j].first != iter->second.pointDataArrays[j].first)  // if the name of the jth field variable is different
            {
              LOG(DEBUG) << "Mesh " << meshName << " cannot be combined with " << vtkPiece_.meshNamesCombinedMeshes << ". Field variable names mismatch for "
                << meshName << " (there is \"" << vtkPiece_.properties.pointDataArrays[j].first << "\" instead of \"" << iter->second.pointDataArrays[j].first << "\")";
              combineMesh = false;
            }
          }
        }

        if (combineMesh)
        {
          VLOG(1) << "Combine mesh " << meshName << " with " << vtkPiece_.meshNamesCombinedMeshes
            << ", add " << iter->second.nPointsLocal << " points, " << iter->second.nCellsLocal << " elements to "
            << vtkPiece_.properties.nPointsLocal << " points, " << vtkPiece_.properties.nCellsLocal << " elements";

          vtkPiece_.properties.nPointsLocal += iter->second.nPointsLocal;
          vtkPiece_.properties.nCellsLocal += iter->second.nCellsLocal;

          vtkPiece_.properties.nPointsGlobal += iter->second.nPointsGlobal;
          vtkPiece_.properties.nCellsGlobal += iter->second.nCellsGlobal;
          vtkPiece_.setVTKValues();
        }
      }
      else
      {
        VLOG(1) << "this is the first " << targetDimensionality << "D mesh";

        // properties are not yet assigned
        vtkPiece_.properties = iter->second; // store properties
        vtkPiece_.setVTKValues();
      }

      VLOG(1) << "combineMesh: " << combineMesh;
      if (combineMesh)
      {
        vtkPiece_.meshNamesCombinedMeshes.insert(meshName);
      }
    }

    LOG(DEBUG) << "vtkPiece_: meshNamesCombinedMeshes: " << vtkPiece_.meshNamesCombinedMeshes << ", properties: " << vtkPiece_.properties
      << ", firstScalarName: " << vtkPiece_.firstScalarName << ", firstVectorName: " << vtkPiece_.firstVectorName;

    Control::PerformanceMeasurement::stop("durationParaview2D3DInit");
  }

  // write combined mesh
  if (!vtkPiece_.meshNamesCombinedMeshes.empty())
  {
    PolyDataPropertiesForMesh &polyDataPropertiesForMesh = vtkPiece_.properties;

    VLOG(1) << "polyDataPropertiesForMesh combined: " << polyDataPropertiesForMesh;

    // if the meshes have the right dimensionality (2D or 3D)
    if (vtkPiece_.properties.dimensionality == targetDimensionality)
    {
      meshNames = vtkPiece_.meshNamesCombinedMeshes;

      std::stringstream filename;
      filename << this->filenameBaseWithNo_ << ".vtu";

      // write actual file
      writeCombinedUnstructuredGridFile<FieldVariablesForOutputWriterType>(
        fieldVariables, vtkPiece_.properties, meshPropertiesUnstructuredGridFile_, vtkPiece_.meshNamesCombinedMeshes, meshPropertiesInitialized, filename.str());
    }
    else
    {
      VLOG(1) << "skip meshes " << vtkPiece_.meshNamesCombinedMeshes << " because " << vtkPiece_.properties.dimensionality << " != " << targetDimensionality;
    }
  }

  // write other single meshes

  // loop over 3D or 2D meshes
  for (std::map<std::string, PolyDataPropertiesForMesh>::iterator meshPropertiesIter = meshPropertiesUnstructuredGridFile_.begin(); meshPropertiesIter != meshPropertiesUnstructuredGridFile_.end(); meshPropertiesIter++)
  {
    PolyDataPropertiesForMesh &polyDataPropertiesForMesh = meshPropertiesIter->second;

    LOG(DEBUG) << "meshNames: " << meshNames << ", next mesh to write: " << meshPropertiesIter->first;

    // if the current mesh was already output within the combined file, skip this mesh
    if (meshNames.find(meshPropertiesIter->first) != meshNames.end())
      continue;

    std::stringstream filename;
    if (meshPropertiesUnstructuredGridFile_.size() > 1)
    {
      filename << this->filenameBaseWithNo_ << "_" << meshPropertiesIter->first << ".vtu";
    }
    else
    {
      std::stringstream filename;
      filename << this->filenameBaseWithNo_ << ".vtu";
    }

    // write actual file
    std::set<std::string> currentMesh;
    currentMesh.insert(meshPropertiesIter->first);
    meshNames.insert(meshPropertiesIter->first);

    writeCombinedUnstructuredGridFile<FieldVariablesForOutputWriterType>(
      fieldVariables, polyDataPropertiesForMesh, meshPropertiesUnstructuredGridFile_, currentMesh, meshPropertiesInitialized, filename.str());

    // for the next meshes, reinitialize, the initialized information can only be shared if there is exactly 1 mesh to write in this method
    meshPropertiesInitialized = false;
  }
}


template<typename FieldVariablesForOutputWriterType>
void Paraview::writeCombinedUnstructuredGridFile(const FieldVariablesForOutputWriterType &fieldVariables, PolyDataPropertiesForMesh &polyDataPropertiesForMesh,
                                                 const std::map<std::string, PolyDataPropertiesForMesh> &meshPropertiesUnstructuredGridFile,
                                                 std::set<std::string> meshNames,
                                                 bool meshPropertiesInitialized, std::string filename)
{
  int targetDimensionality = polyDataPropertiesForMesh.dimensionality;
  bool output3DMeshes = targetDimensionality == 3;

  //! find out name of first scalar and vector field variables
  std::string firstScalarName;
  std::string firstVectorName;

  for (std::vector<std::pair<std::string,int>>::iterator iter = polyDataPropertiesForMesh.pointDataArrays.begin(); iter != polyDataPropertiesForMesh.pointDataArrays.end(); iter++)
  {
    if (iter->second == 1 && firstScalarName.empty())
    {
      firstScalarName = iter->first;
    }
    if (iter->second != 1 && firstVectorName.empty())
    {
      firstVectorName = iter->first;
    }

    if (!firstScalarName.empty() && firstVectorName.empty())
      break;
  }

  // determine filename, broadcast from rank 0
  int filenameLength = filename.length();

  // broadcast length of filename
  MPIUtility::handleReturnValue(MPI_Bcast(&filenameLength, 1, MPI_INT, 0, this->rankSubset_->mpiCommunicator()));

  std::vector<char> receiveBuffer(filenameLength+1, char(0));
  strcpy(receiveBuffer.data(), filename.c_str());
  MPIUtility::handleReturnValue(MPI_Bcast(receiveBuffer.data(), filenameLength, MPI_CHAR, 0, this->rankSubset_->mpiCommunicator()));

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
    nPointsGlobal3D_ = 0;
    MPIUtility::handleReturnValue(MPI_Reduce(&polyDataPropertiesForMesh.nPointsLocal, &nPointsGlobal3D_, 1, MPI_INT, MPI_SUM, 0, this->rankSubset_->mpiCommunicator()), "MPI_Reduce");

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

    VLOG(1) << "nPointsGlobal: " << nPointsGlobal3D_;
  }
  else
  {
    VLOG(1) << "nCellsPreviousRanks3D_, nPointsPreviousRanks3D_, nPointsGlobal3D_ already set to "
      << nCellsPreviousRanks3D_ << ", " << nPointsPreviousRanks3D_ << ", " << nPointsGlobal3D_;
  }

  std::vector<node_no_t> &nNodesLocalWithGhosts = polyDataPropertiesForMesh.nNodesLocalWithGhosts;
  assert(nNodesLocalWithGhosts.size() == targetDimensionality);

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
  }

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

    // loop over meshes
    for (std::string meshName : meshNames)
    {
      int nPointsPreviousRanks3D = connectivityValuesOffset;

      const std::vector<node_no_t> &nNodesLocalWithGhosts = meshPropertiesUnstructuredGridFile.at(meshName).nNodesLocalWithGhosts;

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
      offsetValues[offsetValuesIndexOffset + i] = offsetValuesOffset + (nCellsPreviousRanks3D_[meshName] + i + 1) * nNodesPerCell;    // specifies the end, i.e. one after the last, of the last of nodes for each element
    }

    offsetValuesIndexOffset += meshPropertiesUnstructuredGridFile.at(meshName).nCellsLocal;
    offsetValuesOffset += meshPropertiesUnstructuredGridFile.at(meshName).nCellsLocal * nNodesPerCell;
  }

  VLOG(1) << "offsetValues: " << offsetValues;

  // collect all data for the field variables, organized by field variable names
  std::map<std::string, std::vector<double>> fieldVariableValues;
  ParaviewLoopOverTuple::loopGetNodalValues<FieldVariablesForOutputWriterType>(fieldVariables, meshNames, fieldVariableValues);

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
        pointDataArraysNames << polyDataPropertiesForMesh.pointDataArrays[i].first << " ";
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
    polyDataPropertiesForMesh.pointDataArrays.push_back(std::pair<std::string,int>("partitioning", 1));
  }

  // set data for partitioning field variable
  assert (!fieldVariableValues.empty());
  fieldVariableValues["partitioning"].resize(polyDataPropertiesForMesh.nPointsLocal, (double)this->rankSubset_->ownRankNo());

  // for 2D field variables, add 0 every 2nd entry to generate 3D values, because paraview cannot handle 2D vectors properly
  if (targetDimensionality == 2)
  {
    std::vector<double> buffer;

    // loop over all field variables
    for (std::vector<std::pair<std::string,int>>::iterator pointDataArrayIter = polyDataPropertiesForMesh.pointDataArrays.begin(); pointDataArrayIter != polyDataPropertiesForMesh.pointDataArrays.end(); pointDataArrayIter++)
    {
      // if it is a 2D vector field
      if (pointDataArrayIter->second == 2)
      {

        std::vector<double> &values = fieldVariableValues[pointDataArrayIter->first];

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
  ParaviewLoopOverTuple::loopGetGeometryFieldNodalValues<FieldVariablesForOutputWriterType>(fieldVariables, meshNames, geometryFieldValues);

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

  outputFileParts[outputFilePartNo] << std::string(2, '\t') << "<Piece NumberOfPoints=\"" << nPointsGlobal3D_
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
  for (std::vector<std::pair<std::string,int>>::iterator pointDataArrayIter = polyDataPropertiesForMesh.pointDataArrays.begin(); pointDataArrayIter != polyDataPropertiesForMesh.pointDataArrays.end(); pointDataArrayIter++)
  {
    // paraview cannot handle 2D vector fields without warnings, so set 3D vector fields
    int nComponentsParaview = pointDataArrayIter->second;
    if (nComponentsParaview == 2)
      nComponentsParaview = 3;

    // write normal data element
    outputFileParts[outputFilePartNo] << std::string(4, '\t') << "<DataArray "
        << "Name=\"" << pointDataArrayIter->first << "\" "
        << "type=\"" << (pointDataArrayIter->first == "partitioning"? "Int32" : "Float32") << "\" "
        << "NumberOfComponents=\"" << nComponentsParaview << "\" format=\"" << (binaryOutput_? "binary" : "ascii")
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
  int fieldVariableNo = 0;
  for (std::vector<std::pair<std::string,int>>::iterator pointDataArrayIter = polyDataPropertiesForMesh.pointDataArrays.begin();
      pointDataArrayIter != polyDataPropertiesForMesh.pointDataArrays.end(); pointDataArrayIter++, fieldVariableNo++)
  {
    assert(fieldVariableValues.find(pointDataArrayIter->first) != fieldVariableValues.end());

    VLOG(1) << "write vector for field variable \"" << pointDataArrayIter->first << "\".";

    // write values
    bool writeFloatsAsInt = pointDataArrayIter->first == "partitioning";    // for partitioning, convert float values to integer values for output
    writeCombinedValuesVector(fileHandle, ownRankNo, fieldVariableValues[pointDataArrayIter->first], fieldVariableNo, writeFloatsAsInt);

    // write next xml constructs
    writeAsciiDataShared(fileHandle, ownRankNo, outputFileParts[outputFilePartNo].str());
    outputFilePartNo++;
  }

  VLOG(1) << "write vector for geometry data";

  // write geometry field data
  writeCombinedValuesVector(fileHandle, ownRankNo, geometryFieldValues, fieldVariableNo++);

  // write next xml constructs
  writeAsciiDataShared(fileHandle, ownRankNo, outputFileParts[outputFilePartNo].str());
  outputFilePartNo++;

  // write connectivity values
  writeCombinedValuesVector(fileHandle, ownRankNo, connectivityValues, fieldVariableNo++);

  // write next xml constructs
  writeAsciiDataShared(fileHandle, ownRankNo, outputFileParts[outputFilePartNo].str());
  outputFilePartNo++;

  // write offset values
  writeCombinedValuesVector(fileHandle, ownRankNo, offsetValues, fieldVariableNo++);

  // write next xml constructs
  writeAsciiDataShared(fileHandle, ownRankNo, outputFileParts[outputFilePartNo].str());
  outputFilePartNo++;

  // write types values
  writeCombinedTypesVector(fileHandle, ownRankNo, polyDataPropertiesForMesh.nCellsGlobal, output3DMeshes, fieldVariableNo++);

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
}



} // namespace
