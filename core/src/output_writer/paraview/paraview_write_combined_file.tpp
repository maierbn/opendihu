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

namespace OutputWriter
{

template<typename T>
void Paraview::writeCombinedValuesVector(MPI_File fileHandle, int ownRankNo, const std::vector<T> &values, int identifier)
{
  // fill the write buffer with the local values
  std::string writeBuffer;
  //std::stringstream info;

  if (binaryOutput_)
  {
    VLOG(1) << "Paraview::writeCombinedValuesVector, values: " << values;
    LOG(DEBUG) << "rankSubset: " << *this->rankSubset_;

    int localValuesSize = values.size() * sizeof(float);  // number of bytes
    int nLocalValues = values.size() + (ownRankNo == 0? 1 : 0);

    // initialize cached values
    if (globalValuesSize_.size() < identifier+1)
    {
      globalValuesSize_.resize(identifier+1);

      // gather data length of total vector to rank 0
      int globalValuesSize = 0;
      MPIUtility::handleReturnValue(MPI_Reduce(&localValuesSize, &globalValuesSize, 1, MPI_INT,
                                                MPI_SUM, 0, this->rankSubset_->mpiCommunicator()));

      globalValuesSize_[identifier] = globalValuesSize;  // value only set on rank 0, all other ranks have value 0

      // determine number of previous values
      nPreviousValues_.resize(identifier+1);
      int nPreviousValues = 0;
      MPIUtility::handleReturnValue(MPI_Exscan(&nLocalValues, &nPreviousValues, 1, MPI_INT, MPI_SUM, this->rankSubset_->mpiCommunicator()), "MPI_Exscan");
      nPreviousValues_[identifier] = nPreviousValues;

      VLOG(1) << "reduce data length of total vector, localValuesSize: " << localValuesSize << ", global size: " << globalValuesSize << ", nPreviousValues: " << nPreviousValues;
    }

    std::list<int32_t> valuesVector;

    // convert float values to int32 and insert into valuesVector
    if (std::is_same<T,double>::value || std::is_same<T,float>::value)
    {
      for (int i = 0; i < values.size(); i++)
      {
        union
        {
          float f;
          int32_t j;
        };
        f = values[i];
        valuesVector.push_back(j);
      }
    }
    else
    {
      // copy all int values
      std::copy(values.begin(), values.end(), std::back_inserter(valuesVector));
    }

    // on rank 0 prepend the encoded total length of the data vector in bytes
    if (ownRankNo == 0)
    {
      valuesVector.push_front(globalValuesSize_[identifier]);
    }

    // send first value to previous rank
    if (valuesVector.size() < 1)
    {
      LOG(FATAL) << "No values for rank " << ownRankNo << " in paraview binary output. This is not supported for binary output, so set binary to False.";
    }
    float firstValue = (float)valuesVector.front();
    MPI_Request sendRequest;
    MPI_Request receiveRequest;
    if (ownRankNo != 0)
    {
      MPIUtility::handleReturnValue(MPI_Isend(&firstValue, 1, MPI_FLOAT, ownRankNo-1, 0, this->rankSubset_->mpiCommunicator(), &sendRequest), "MPI_Isend");
    }

    // receive the next ranks first value
    float firstValueNextRank = 0;
    if (ownRankNo != this->rankSubset_->size()-1)
    {
      MPIUtility::handleReturnValue(MPI_Irecv(&firstValueNextRank, 1, MPI_FLOAT, ownRankNo+1, 0, this->rankSubset_->mpiCommunicator(), &receiveRequest), "MPI_Irecv");
      MPI_Status status;
      MPIUtility::handleReturnValue(MPI_Wait(&receiveRequest, &status), "MPI_Wait", &status);
    }
    if (ownRankNo != 0)
    {
      MPI_Status status;
      MPIUtility::handleReturnValue(MPI_Wait(&sendRequest, &status), "MPI_Wait", &status);
    }

    VLOG(1) << "firstValueNextRank: " << firstValueNextRank << ", nPreviousValues_[" << identifier << "]: " << nPreviousValues_[identifier] << ", %3=" << (nPreviousValues_[identifier]%3);

    // the end of the own values are contained in the own written data plus some part of the first value of the next rank
    if (nPreviousValues_[identifier] % 3 == 0)
    {
      valuesVector.push_back(firstValueNextRank);
      writeBuffer = Paraview::encodeBase64Int(valuesVector.begin(), valuesVector.end(), false);  //without leading dataset size

      // remove last bytes that encode the additional last value
      int offset = 0;
      int nBits = (std::distance(valuesVector.begin(), valuesVector.end())-1)*4*8 - offset;
      int nBytes = nBits / 6;
      if (nBits % 6 != 0)
      {
        nBytes = int(nBits / 6 + 1);
      }
      writeBuffer = writeBuffer.substr(0, nBytes);

      VLOG(1) << "offset = 0 bits";
      //info << "offset = 0 bits";

    }
    else if (nPreviousValues_[identifier] % 3 == 1)
    {
      // offset = 4 bits (4 bits were written by the previous rank)
      // add dummy value (4 bytes)
      valuesVector.push_front(0);
      valuesVector.push_back(firstValueNextRank);
      std::string stringData = Paraview::encodeBase64Int(valuesVector.begin(), valuesVector.end(), false);  //without leading dataset size

      VLOG(1) << "offset = 4 bits [" << stringData << "]";
      //info << "offset = 4 bits [" << stringData << "]";

      // remove last bytes that encode the additional last value
      int offset = 4;
      int nBits = (std::distance(valuesVector.begin(), valuesVector.end())-2)*4*8 - offset;
      int nBytes = nBits / 6;
      if (nBits % 6 != 0)
      {
        nBytes = int(nBits / 6 + 1);
      }

      stringData = stringData.substr(6);    // discard first six bytes (discards padding 4 bytes and gives offset of 4 bits (6*6=36 = 4*8+4))
      if (ownRankNo != this->rankSubset_->size()-1)
      {
        stringData = stringData.substr(0, nBytes);
      }

      VLOG(1) << "padding removed [" << stringData << "] (size=" << this->rankSubset_->size() << ")";
      //info << ", padding removed [" << stringData << "] nBits: " << (std::distance(valuesVector.begin(), valuesVector.end())-2)
      //  << "*4*8 - " << offset << " = " << nBits << ", nBits % 6 = " << nBits % 6 << ", nBytes = " << nBytes << "  ";

      writeBuffer = stringData;
    }
    else if (nPreviousValues_[identifier] % 3 == 2)
    {
      // offset = 2 bits (2 bits were written by the previous rank)
      // add two dummy values (8 bytes)
      valuesVector.push_front(0);
      valuesVector.push_front(0);
      valuesVector.push_back(firstValueNextRank);
      std::string stringData = Paraview::encodeBase64Int(valuesVector.begin(), valuesVector.end(), false);  //without leading dataset size

      VLOG(1) << "offset = 2 bits [" << stringData << "]";
      //info << "offset = 2 bits [" << stringData << "]";

      // remove last bytes that encode the additional last value
      int offset = 2;
      int nBits = (std::distance(valuesVector.begin(), valuesVector.end())-3)*4*8 - offset;
      int nBytes = nBits / 6;

      if (nBits % 6 != 0)
      {
        nBytes = int(nBits / 6 + 1);
      }

      stringData = stringData.substr(11);  // discard first 11 bytes (discards padding 2*4=8 bytes and gives offset of 2 bits (11*6=66 = 8*8+2))

      if (ownRankNo != this->rankSubset_->size()-1)
      {
        stringData = stringData.substr(0, nBytes);
      }

      VLOG(1) << "padding removed [" << stringData << "]";
      //info << ", padding removed [" << stringData << "]";

      writeBuffer = stringData;
    }

    if (VLOG_IS_ON(1))
    {
      VLOG(1) << " rank " << ownRankNo << ", values: ";
      for (std::list<int32_t>::iterator iter=valuesVector.begin(); iter != valuesVector.end(); iter++)
        VLOG(1) << *iter << " ";
    }

    // info << " (rank " << ownRankNo << ", values:";
    // for (std::list<int32_t>::iterator iter=valuesVector.begin(); iter != valuesVector.end(); iter++)
    //   info << " " << *iter;
    // info << ") ";

    // determine global number of bytes and make sure that it is a multiple of 4 by adding padding "="s
    int nBytesLocal = writeBuffer.length();
    int nBytesGlobal = 0;
    MPIUtility::handleReturnValue(MPI_Reduce(&nBytesLocal, &nBytesGlobal, 1, MPI_INT,
                                             MPI_SUM, this->rankSubset_->size()-1, this->rankSubset_->mpiCommunicator()));

    // add base64 padding with "=" characters such that global size is a factor of 4
    if (ownRankNo == this->rankSubset_->size()-1)
    {
      // for the last rank remove base64 padding with '='s
      Base64::StripPadding(&writeBuffer);

      if (nBytesGlobal % 4 == 1)
      {
        // this should not happen
        writeBuffer += std::string("===");

        LOG(FATAL) << "Base64 encoding in parallel has a bug, padding needs 3 '=' signs, but only 1 or 2 should be ever needed";
      }
      else if (nBytesGlobal % 4 == 2)
      {
        writeBuffer += std::string("==");
      }
      else if (nBytesGlobal % 4 == 3)
      {
        writeBuffer += std::string("=");
      }
    }
  }
  else
  {
    writeBuffer = Paraview::convertToAscii(values, fixedFormat_);
  }

  VLOG(1) << "writeBuffer [" << writeBuffer << "]";

  // writeBuffer = info.str() + writeBuffer;

  // add line break after each process (for debugging)
  if (ownRankNo != this->rankSubset_->size()-1 && false)
  {
    writeBuffer += std::string("\n");
    writeBuffer += std::string(5,'\t');
  }

  MPIUtility::handleReturnValue(MPI_File_write_ordered(fileHandle, writeBuffer.c_str(), writeBuffer.length(), MPI_BYTE, MPI_STATUS_IGNORE), "MPI_File_write_ordered");
}

template<typename OutputFieldVariablesType>
void Paraview::writePolyDataFile(const OutputFieldVariablesType &fieldVariables, std::set<std::string> &combinedMeshesOut)
{
  // output a *.vtp file which contains 1D meshes, if there are any

  // collect the size data that is needed to compute offsets for parallel file output
  std::map<std::string, PolyDataPropertiesForMesh> meshProperties;
  ParaviewLoopOverTuple::loopCollectMeshProperties<OutputFieldVariablesType>(fieldVariables, meshProperties);

  VLOG(1) << "writePolyDataFile on rankSubset_: " << *this->rankSubset_;
  assert(this->rankSubset_);

  VLOG(1) << "meshProperties: " << meshProperties;

  /*
   PolyDataPropertiesForMesh:
  int dimensionality;    ///< D=1: object is a VTK "Line" and can be written to a combined vtp file
  global_no_t nPoints;   ///< the number of points needed for representing the mesh, note that some points may be
                         ///< duplicated because they are ghosts on one process and non-ghost on the other. This is intended.
  global_no_t nCells;    ///< the number of VTK "cells", i.e. "Lines", which is the opendihu number of "1D elements"

  std::vector<std::pair<std::string,int>> pointDataArrays;   ///< <name,nComponents> of PointData DataArray elements
  */

  /* one VTKPiece is the XML element that will be output as <Piece></Piece>. It is created from one or multiple opendihu meshes
   */
  struct VTKPiece
  {
    std::set<std::string> meshNamesCombinedMeshes;   ///< the meshNames of the combined meshes, or only one meshName if it is not a merged mesh
    PolyDataPropertiesForMesh properties;   ///< the properties of the merged mesh

    std::string firstScalarName;   ///< name of the first scalar field variable of the mesh
    std::string firstVectorName;   ///< name of the first non-scalar field variable of the mesh

    //! constructor, initialize nPoints and nCells to 0
    VTKPiece()
    {
      properties.nPointsLocal = 0;
      properties.nCellsLocal = 0;
      properties.nPointsGlobal = 0;
      properties.nCellsGlobal = 0;
    }

    //! assign the correct values to firstScalarName and firstVectorName, only if properties has been set
    void setVTKValues()
    {
      // set values for firstScalarName and firstVectorName from the values in pointDataArrays
      for (auto pointDataArray : properties.pointDataArrays)
      {
        if (firstScalarName == "" && pointDataArray.second == 1)
          firstScalarName = pointDataArray.first;
        if (firstVectorName == "" && pointDataArray.second != 1)
          firstVectorName = pointDataArray.first;
      }
    }
  } vtkPiece;

  // parse the collected properties of the meshes that will be output to the file
  for (std::map<std::string,PolyDataPropertiesForMesh>::iterator iter = meshProperties.begin(); iter != meshProperties.end(); iter++)
  {
    std::string meshName = iter->first;

    // do not combine meshes other than 1D meshes
    if (iter->second.dimensionality != 1)
      continue;

    // check if this mesh should be combined with other meshes
    bool combineMesh = true;

    // check if mesh can be merged into previous meshes
    if (!vtkPiece.properties.pointDataArrays.empty())   // if properties are already assigned by an earlier mesh
    {
      if (vtkPiece.properties.pointDataArrays.size() != iter->second.pointDataArrays.size())
      {
        LOG(DEBUG) << "Mesh " << meshName << " cannot be combined with " << vtkPiece.meshNamesCombinedMeshes << ". Number of field variables mismatches for "
          << meshName << " (is " << iter->second.pointDataArrays.size() << " instead of " << vtkPiece.properties.pointDataArrays.size() << ")";
        combineMesh = false;
      }
      else
      {
        for (int j = 0; j < iter->second.pointDataArrays.size(); j++)
        {
          if (vtkPiece.properties.pointDataArrays[j].first != iter->second.pointDataArrays[j].first)  // if the name of the jth field variable is different
          {
            LOG(DEBUG) << "Mesh " << meshName << " cannot be combined with " << vtkPiece.meshNamesCombinedMeshes << ". Field variable names mismatch for "
              << meshName << " (there is \"" << vtkPiece.properties.pointDataArrays[j].first << "\" instead of \"" << iter->second.pointDataArrays[j].first << "\")";
            combineMesh = false;
          }
        }
      }

      if (combineMesh)
      {
        LOG(DEBUG) << "Combine mesh " << meshName << " with " << vtkPiece.meshNamesCombinedMeshes
          << ", add " << iter->second.nPointsLocal << " points, " << iter->second.nCellsLocal << " elements to "
          << vtkPiece.properties.nPointsLocal << " points, " << vtkPiece.properties.nCellsLocal << " elements";

        vtkPiece.properties.nPointsLocal += iter->second.nPointsLocal;
        vtkPiece.properties.nCellsLocal += iter->second.nCellsLocal;

        vtkPiece.properties.nPointsGlobal += iter->second.nPointsGlobal;
        vtkPiece.properties.nCellsGlobal += iter->second.nCellsGlobal;
        vtkPiece.setVTKValues();
      }
    }
    else
    {
      VLOG(1) << "this is the first 1D mesh";

      // properties are not yet assigned
      vtkPiece.properties = iter->second; // store properties
      vtkPiece.setVTKValues();
    }

    VLOG(1) << "combineMesh: " << combineMesh;
    if (combineMesh)
    {
      vtkPiece.meshNamesCombinedMeshes.insert(meshName);
    }
  }

  LOG(DEBUG) << "vtkPiece: meshNamesCombinedMeshes: " << vtkPiece.meshNamesCombinedMeshes << ", properties: " << vtkPiece.properties
    << ", firstScalarName: " << vtkPiece.firstScalarName << ", firstVectorName: " << vtkPiece.firstVectorName;

  combinedMeshesOut = vtkPiece.meshNamesCombinedMeshes;

  // add field variable "partitioning" with 1 component
  vtkPiece.properties.pointDataArrays.push_back(std::pair<std::string,int>("partitioning", 1));

  // determine filename, broadcast from rank 0
  std::stringstream filename;
  filename << this->filenameBaseWithNo_ << ".vtp";
  int filenameLength = filename.str().length();

  // broadcast length of filename
  MPIUtility::handleReturnValue(MPI_Bcast(&filenameLength, 1, MPI_INT, 0, this->rankSubset_->mpiCommunicator()));

  std::vector<char> receiveBuffer(filenameLength+1, char(0));
  strcpy(receiveBuffer.data(), filename.str().c_str());
  MPIUtility::handleReturnValue(MPI_Bcast(receiveBuffer.data(), filenameLength, MPI_CHAR, 0, this->rankSubset_->mpiCommunicator()));

  std::string filenameStr(receiveBuffer.begin(), receiveBuffer.end());

  // remove file if it exists, synchronization afterwards by MPI calls, that is why the remove call is already here
  assert(this->rankSubset_);
  int ownRankNo = this->rankSubset_->ownRankNo();
  if (ownRankNo == 0)
  {
    // open file to ensure that directory exists and file is writable
    std::ofstream file = Generic::openFile(filenameStr);

    // close and delete file
    file.close();
    std::remove(filenameStr.c_str());
  }

  // exchange information about offset in terms of nCells and nPoints
  int nCellsPreviousRanks = 0;
  int nPointsPreviousRanks = 0;
  int nPointsGlobal = 0;
  int nLinesGlobal = 0;

  MPIUtility::handleReturnValue(MPI_Exscan(&vtkPiece.properties.nCellsLocal, &nCellsPreviousRanks, 1, MPI_INT, MPI_SUM, this->rankSubset_->mpiCommunicator()), "MPI_Exscan");
  MPIUtility::handleReturnValue(MPI_Exscan(&vtkPiece.properties.nPointsLocal, &nPointsPreviousRanks, 1, MPI_INT, MPI_SUM, this->rankSubset_->mpiCommunicator()), "MPI_Exscan");
  MPIUtility::handleReturnValue(MPI_Reduce(&vtkPiece.properties.nPointsLocal, &nPointsGlobal, 1, MPI_INT, MPI_SUM, 0, this->rankSubset_->mpiCommunicator()), "MPI_Reduce");
  MPIUtility::handleReturnValue(MPI_Reduce(&vtkPiece.properties.nCellsLocal, &nLinesGlobal, 1, MPI_INT, MPI_SUM, 0, this->rankSubset_->mpiCommunicator()), "MPI_Reduce");

  // get local data values
  // setup connectivity array
  std::vector<int> connectivityValues(2*vtkPiece.properties.nCellsLocal);
  for (int i = 0; i < vtkPiece.properties.nCellsLocal; i++)
  {
    connectivityValues[2*i + 0] = nPointsPreviousRanks + i;
    connectivityValues[2*i + 1] = nPointsPreviousRanks + i+1;
  }

  // setup offset array
  std::vector<int> offsetValues(vtkPiece.properties.nCellsLocal);
  for (int i = 0; i < vtkPiece.properties.nCellsLocal; i++)
  {
    offsetValues[i] = 2*nCellsPreviousRanks + 2*i + 1;
  }

  // collect all data for the field variables, organized by field variable names
  std::map<std::string, std::vector<double>> fieldVariableValues;
  ParaviewLoopOverTuple::loopGetNodalValues<OutputFieldVariablesType>(fieldVariables, vtkPiece.meshNamesCombinedMeshes, fieldVariableValues);

  assert (!fieldVariableValues.empty());
  fieldVariableValues["partitioning"].resize(vtkPiece.properties.nPointsLocal, (double)this->rankSubset_->ownRankNo());

  // if next assertion will fail, output why for debugging
  if (fieldVariableValues.size() != vtkPiece.properties.pointDataArrays.size())
  {
    LOG(DEBUG) << "n field variable values: " << fieldVariableValues.size() << ", n point data arrays: "
      << vtkPiece.properties.pointDataArrays.size();
    LOG(DEBUG) << "vtkPiece.meshNamesCombinedMeshes: " << vtkPiece.meshNamesCombinedMeshes;
    std::stringstream pointDataArraysNames;
    for (int i = 0; i < vtkPiece.properties.pointDataArrays.size(); i++)
    {
      pointDataArraysNames << vtkPiece.properties.pointDataArrays[i].first << " ";
    }
    LOG(DEBUG) << "pointDataArraysNames: " <<  pointDataArraysNames.str();
  }

  assert(fieldVariableValues.size() == vtkPiece.properties.pointDataArrays.size());

  // collect all data for the geometry field variable
  std::vector<double> geometryFieldValues;
  ParaviewLoopOverTuple::loopGetGeometryFieldNodalValues<OutputFieldVariablesType>(fieldVariables, vtkPiece.meshNamesCombinedMeshes, geometryFieldValues);

  // only continue if there is data to reduce
  if (vtkPiece.meshNamesCombinedMeshes.empty())
  {
    LOG(ERROR) << "There are no 1D meshes that could be combined, but Paraview output with combineFiles=True was specified. \n(This only works for 1D meshes.)";
  }

  LOG(DEBUG) << "Combined mesh from " << vtkPiece.meshNamesCombinedMeshes;

  int nOutputFileParts = 4 + vtkPiece.properties.pointDataArrays.size();

  // create the basic structure of the output file
  std::vector<std::stringstream> outputFileParts(nOutputFileParts);
  int outputFilePartNo = 0;
  outputFileParts[outputFilePartNo] << "<?xml version=\"1.0\"?>" << std::endl
    << "<!-- " << DihuContext::versionText() << " " << DihuContext::metaText() << "-->" << std::endl
    << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\">" << std::endl    // intel cpus are LittleEndian
    << std::string(1, '\t') << "<PolyData>" << std::endl;

  outputFileParts[outputFilePartNo] << std::string(2, '\t') << "<Piece NumberOfPoints=\"" << nPointsGlobal << "\" NumberOfVerts=\"0\" "
    << "NumberOfLines=\"" << nLinesGlobal << "\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">" << std::endl
    << std::string(3, '\t') << "<PointData";

  if (vtkPiece.firstScalarName != "")
  {
    outputFileParts[outputFilePartNo] << " Scalars=\"" << vtkPiece.firstScalarName << "\"";
  }
  if (vtkPiece.firstVectorName != "")
  {
    outputFileParts[outputFilePartNo] << " Vectors=\"" << vtkPiece.firstVectorName << "\"";
  }
  outputFileParts[outputFilePartNo] << ">" << std::endl;

  // loop over field variables (PointData)
  for (std::vector<std::pair<std::string,int>>::iterator pointDataArrayIter = vtkPiece.properties.pointDataArrays.begin(); pointDataArrayIter != vtkPiece.properties.pointDataArrays.end(); pointDataArrayIter++)
  {
    // write normal data element
    outputFileParts[outputFilePartNo] << std::string(4, '\t') << "<DataArray "
        << "Name=\"" << pointDataArrayIter->first << "\" "
        << "type=\"Float32\" "
        << "NumberOfComponents=\"" << pointDataArrayIter->second << "\" format=\"" << (binaryOutput_? "binary" : "ascii")
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
    << std::string(3, '\t') << "<Verts></Verts>" << std::endl
    << std::string(3, '\t') << "<Lines>" << std::endl
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
    << std::string(3, '\t') << "</Lines>" << std::endl
    << std::string(3, '\t') << "<Strips></Strips>" << std::endl
    << std::string(3, '\t') << "<Polys></Polys>" << std::endl
    << std::string(2, '\t') << "</Piece>" << std::endl
    << std::string(1, '\t') << "</PolyData>" << std::endl
    << "</VTKFile>" << std::endl;

  assert(outputFilePartNo+1 == nOutputFileParts);

  // loop over output file parts and collect the missing data for the own rank

  VLOG(1) << "outputFileParts:";

  for (std::vector<std::stringstream>::iterator iter = outputFileParts.begin(); iter != outputFileParts.end(); iter++)
  {
    VLOG(1) << "  " << iter->str();
  }

  LOG(DEBUG) << "open MPI file \"" << filenameStr << "\".";

  // open file
  MPI_File fileHandle;
  MPIUtility::handleReturnValue(MPI_File_open(this->rankSubset_->mpiCommunicator(), filenameStr.c_str(),
                                              //MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_UNIQUE_OPEN,
                                              MPI_MODE_WRONLY | MPI_MODE_CREATE,
                                              MPI_INFO_NULL, &fileHandle), "MPI_File_open");

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
  for (std::vector<std::pair<std::string,int>>::iterator pointDataArrayIter = vtkPiece.properties.pointDataArrays.begin();
       pointDataArrayIter != vtkPiece.properties.pointDataArrays.end(); pointDataArrayIter++, fieldVariableNo++)
  {
    assert(fieldVariableValues.find(pointDataArrayIter->first) != fieldVariableValues.end());

    // write values
    writeCombinedValuesVector(fileHandle, ownRankNo, fieldVariableValues[pointDataArrayIter->first], fieldVariableNo);

    // write next xml constructs
    writeAsciiDataShared(fileHandle, ownRankNo, outputFileParts[outputFilePartNo].str());
    outputFilePartNo++;
  }

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

  MPIUtility::handleReturnValue(MPI_File_close(&fileHandle), "MPI_File_close");
}

} // namespace
