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
#include "control/performance_measurement.h"

namespace OutputWriter
{

template<typename T>
void Paraview::writeCombinedValuesVector(MPI_File fileHandle, int ownRankNo, const std::vector<T> &values, int identifier, bool writeFloatsAsInt)
{
  // fill the write buffer with the local values
  std::string writeBuffer;
  //std::stringstream info;

  if (binaryOutput_)
  {
    VLOG(1) << "Paraview::writeCombinedValuesVector, " << values.size() << " values: " << values;
    VLOG(1) << "rankSubset: " << *this->rankSubset_;

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
      if (writeFloatsAsInt)
      {
        // values are float type but should be converted to integer values
        for (int i = 0; i < values.size(); i++)
        {
          int32_t integerValue = (int32_t)(round(values[i]));
          valuesVector.push_back(integerValue);
        }
      }
      else
      {
        // values are float type, put 4 byte data per float in list of int32_t datatype
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
    }
    else
    {
      // values are already int32_t
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
      writeBuffer = Paraview::encodeBase64Int32(valuesVector.begin(), valuesVector.end(), false);  //without leading dataset size

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
      std::string stringData = Paraview::encodeBase64Int32(valuesVector.begin(), valuesVector.end(), false);  //without leading dataset size

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
      std::string stringData = Paraview::encodeBase64Int32(valuesVector.begin(), valuesVector.end(), false);  //without leading dataset size

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

    if (VLOG_IS_ON(3))
    {
      VLOG(3) << " rank " << ownRankNo << ", values: ";
      for (std::list<int32_t>::iterator iter=valuesVector.begin(); iter != valuesVector.end(); iter++)
        VLOG(3) << *iter << " ";
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

  VLOG(2) << "writeBuffer [" << writeBuffer << "]";

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
void Paraview::writePolyDataFile(const OutputFieldVariablesType &fieldVariables, std::set<std::string> &meshNames)
{
  // output a *.vtp file which contains 1D meshes, if there are any

  bool meshPropertiesInitialized = !meshPropertiesPolyDataFile_.empty();

  if (!meshPropertiesInitialized)
  {
    Control::PerformanceMeasurement::start("durationParaview1DInit");

    // collect the size data that is needed to compute offsets for parallel file output
    ParaviewLoopOverTuple::loopCollectMeshProperties<OutputFieldVariablesType>(fieldVariables, meshPropertiesPolyDataFile_);

    Control::PerformanceMeasurement::stop("durationParaview1DInit");
  }

  VLOG(1) << "writePolyDataFile on rankSubset_: " << *this->rankSubset_;
  assert(this->rankSubset_);

  VLOG(1) << "meshPropertiesPolyDataFile_: " << meshPropertiesPolyDataFile_;
  /*
  struct PolyDataPropertiesForMesh
  {
    int dimensionality;    ///< D=1: object is a VTK "Line", D=2, D=3: object is a VTK "Poly"
    global_no_t nPointsLocal;   ///< the number of points needed for representing the mesh, local value of rank
    global_no_t nCellsLocal;    ///< the number of VTK "cells", i.e. "Lines" or "Polys", which is the opendihu number of "elements", local value of rank
    global_no_t nPointsGlobal;   ///< the number of points needed for representing the mesh, global value of all rank
    global_no_t nCellsGlobal;    ///< the number of VTK "cells", i.e. "Lines" or "Polys", which is the opendihu number of "elements", global value of all ranks

    std::vector<std::pair<std::string,int>> pointDataArrays;   ///< <name,nComponents> of PointData DataArray elements
  };
*/

  /* one VTKPiece is the XML element that will be output as <Piece></Piece>. It is created from one or multiple opendihu meshes
   */
  /*
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
      properties.dimensionality = 0;
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
  */

  if (!meshPropertiesInitialized)
  {
    Control::PerformanceMeasurement::start("durationParaview1DInit");

    // parse the collected properties of the meshes that will be output to the file
    for (std::map<std::string,PolyDataPropertiesForMesh>::iterator iter = meshPropertiesPolyDataFile_.begin(); iter != meshPropertiesPolyDataFile_.end(); iter++)
    {
      std::string meshName = iter->first;

      // do not combine meshes other than 1D meshes
      if (iter->second.dimensionality != 1)
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
        VLOG(1) << "this is the first 1D mesh";

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

    Control::PerformanceMeasurement::stop("durationParaview1DInit");
  }

  meshNames = vtkPiece_.meshNamesCombinedMeshes;

  // if there are no 1D meshes, return
  if (meshNames.empty())
    return;

  if (!meshPropertiesInitialized)
  {
    // add field variable "partitioning" with 1 component
    vtkPiece_.properties.pointDataArrays.push_back(std::pair<std::string,int>("partitioning", 1));
  }

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
    std::ofstream file;
    Generic::openFile(file, filenameStr);

    // close and delete file
    file.close();
    std::remove(filenameStr.c_str());
  }

  if (!meshPropertiesInitialized)
  {
    // exchange information about offset in terms of nCells and nPoints
    nCellsPreviousRanks1D_ = 0;
    nPointsPreviousRanks1D_ = 0;
    nPointsGlobal1D_ = 0;
    nLinesGlobal1D_ = 0;

    Control::PerformanceMeasurement::start("durationParaview1DInit");
    Control::PerformanceMeasurement::start("durationParaview1DReduction");
    MPIUtility::handleReturnValue(MPI_Exscan(&vtkPiece_.properties.nCellsLocal, &nCellsPreviousRanks1D_, 1, MPI_INT, MPI_SUM, this->rankSubset_->mpiCommunicator()), "MPI_Exscan");
    MPIUtility::handleReturnValue(MPI_Exscan(&vtkPiece_.properties.nPointsLocal, &nPointsPreviousRanks1D_, 1, MPI_INT, MPI_SUM, this->rankSubset_->mpiCommunicator()), "MPI_Exscan");
    MPIUtility::handleReturnValue(MPI_Reduce(&vtkPiece_.properties.nPointsLocal, &nPointsGlobal1D_, 1, MPI_INT, MPI_SUM, 0, this->rankSubset_->mpiCommunicator()), "MPI_Reduce");
    MPIUtility::handleReturnValue(MPI_Reduce(&vtkPiece_.properties.nCellsLocal, &nLinesGlobal1D_, 1, MPI_INT, MPI_SUM, 0, this->rankSubset_->mpiCommunicator()), "MPI_Reduce");
    Control::PerformanceMeasurement::stop("durationParaview1DReduction");
    Control::PerformanceMeasurement::stop("durationParaview1DInit");
  }

  // get local data values
  // setup connectivity array
  std::vector<int> connectivityValues(2*vtkPiece_.properties.nCellsLocal);
  for (int i = 0; i < vtkPiece_.properties.nCellsLocal; i++)
  {
    connectivityValues[2*i + 0] = nPointsPreviousRanks1D_ + i;
    connectivityValues[2*i + 1] = nPointsPreviousRanks1D_ + i+1;
  }

  // setup offset array
  std::vector<int> offsetValues(vtkPiece_.properties.nCellsLocal);
  for (int i = 0; i < vtkPiece_.properties.nCellsLocal; i++)
  {
    offsetValues[i] = 2*nCellsPreviousRanks1D_ + 2*i + 1;
  }

  // collect all data for the field variables, organized by field variable names
  std::map<std::string, std::vector<double>> fieldVariableValues;
  ParaviewLoopOverTuple::loopGetNodalValues<OutputFieldVariablesType>(fieldVariables, vtkPiece_.meshNamesCombinedMeshes, fieldVariableValues);

  assert (!fieldVariableValues.empty());
  fieldVariableValues["partitioning"].resize(vtkPiece_.properties.nPointsLocal, (double)this->rankSubset_->ownRankNo());

  // if next assertion will fail, output why for debugging
  if (fieldVariableValues.size() != vtkPiece_.properties.pointDataArrays.size())
  {
    LOG(DEBUG) << "n field variable values: " << fieldVariableValues.size() << ", n point data arrays: "
      << vtkPiece_.properties.pointDataArrays.size();
    LOG(DEBUG) << "vtkPiece_.meshNamesCombinedMeshes: " << vtkPiece_.meshNamesCombinedMeshes;
    std::stringstream pointDataArraysNames;
    for (int i = 0; i < vtkPiece_.properties.pointDataArrays.size(); i++)
    {
      pointDataArraysNames << vtkPiece_.properties.pointDataArrays[i].first << " ";
    }
    LOG(DEBUG) << "pointDataArraysNames: " <<  pointDataArraysNames.str();
  }

  assert(fieldVariableValues.size() == vtkPiece_.properties.pointDataArrays.size());

  // collect all data for the geometry field variable
  std::vector<double> geometryFieldValues;
  ParaviewLoopOverTuple::loopGetGeometryFieldNodalValues<OutputFieldVariablesType>(fieldVariables, vtkPiece_.meshNamesCombinedMeshes, geometryFieldValues);

  // only continue if there is data to reduce
  if (vtkPiece_.meshNamesCombinedMeshes.empty())
  {
    LOG(ERROR) << "There are no 1D meshes that could be combined, but Paraview output with combineFiles=True was specified. \n(This only works for 1D meshes.)";
  }

  LOG(DEBUG) << "Combined mesh from " << vtkPiece_.meshNamesCombinedMeshes;

  int nOutputFileParts = 4 + vtkPiece_.properties.pointDataArrays.size();

  // create the basic structure of the output file
  std::vector<std::stringstream> outputFileParts(nOutputFileParts);
  int outputFilePartNo = 0;
  outputFileParts[outputFilePartNo] << "<?xml version=\"1.0\"?>" << std::endl
    << "<!-- " << DihuContext::versionText() << " " << DihuContext::metaText()
    << ", currentTime: " << this->currentTime_ << ", timeStepNo: " << this->timeStepNo_ << " -->" << std::endl
    << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\">" << std::endl    // intel cpus are LittleEndian
    << std::string(1, '\t') << "<PolyData>" << std::endl;

  outputFileParts[outputFilePartNo] << std::string(2, '\t') << "<Piece NumberOfPoints=\"" << nPointsGlobal1D_ << "\" NumberOfVerts=\"0\" "
    << "NumberOfLines=\"" << nLinesGlobal1D_ << "\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">" << std::endl
    << std::string(3, '\t') << "<PointData";

  if (vtkPiece_.firstScalarName != "")
  {
    outputFileParts[outputFilePartNo] << " Scalars=\"" << vtkPiece_.firstScalarName << "\"";
  }
  if (vtkPiece_.firstVectorName != "")
  {
    outputFileParts[outputFilePartNo] << " Vectors=\"" << vtkPiece_.firstVectorName << "\"";
  }
  outputFileParts[outputFilePartNo] << ">" << std::endl;

  // loop over field variables (PointData)
  for (std::vector<std::pair<std::string,int>>::iterator pointDataArrayIter = vtkPiece_.properties.pointDataArrays.begin(); pointDataArrayIter != vtkPiece_.properties.pointDataArrays.end(); pointDataArrayIter++)
  {
    // write normal data element
    outputFileParts[outputFilePartNo] << std::string(4, '\t') << "<DataArray "
        << "Name=\"" << pointDataArrayIter->first << "\" "
        << "type=\"" << (pointDataArrayIter->first == "partitioning"? "Int32" : "Float32") << "\" "
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

  Control::PerformanceMeasurement::start("durationParaview1DWrite");

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
  for (std::vector<std::pair<std::string,int>>::iterator pointDataArrayIter = vtkPiece_.properties.pointDataArrays.begin();
       pointDataArrayIter != vtkPiece_.properties.pointDataArrays.end(); pointDataArrayIter++, fieldVariableNo++)
  {
    assert(fieldVariableValues.find(pointDataArrayIter->first) != fieldVariableValues.end());

    // write values
    bool writeFloatsAsInt = pointDataArrayIter->first == "partitioning";    // for partitioning, convert float values to integer values for output
    writeCombinedValuesVector(fileHandle, ownRankNo, fieldVariableValues[pointDataArrayIter->first], fieldVariableNo, writeFloatsAsInt);

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

  Control::PerformanceMeasurement::stop("durationParaview1DWrite");

  MPIUtility::handleReturnValue(MPI_File_close(&fileHandle), "MPI_File_close");
}

template<typename OutputFieldVariablesType>
void Paraview::writeCombinedUnstructuredGridFile(const OutputFieldVariablesType &fieldVariables, std::set<std::string> &meshNames,
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
    ParaviewLoopOverTuple::loopCollectMeshProperties<OutputFieldVariablesType>(fieldVariables, meshPropertiesUnstructuredGridFile_);

    Control::PerformanceMeasurement::stop("durationParaview3DInit");
  }
  else LOG(DEBUG) << "meshPropertiesUnstructuredGridFile_ already initialized";

  VLOG(1) << "writeCombinedUnstructuredGridFile on rankSubset_: " << *this->rankSubset_;
  assert(this->rankSubset_);

  VLOG(1) << "meshPropertiesUnstructuredGridFile_: " << meshPropertiesUnstructuredGridFile_;
  /*
  struct PolyDataPropertiesForMesh
  {
    int dimensionality;    ///< D=1: object is a VTK "Line", D=2, D=3: object should be represented by an unstructured grid
    global_no_t nPointsLocal;   ///< the number of points needed for representing the mesh, local value of rank
    global_no_t nCellsLocal;    ///< the number of VTK "cells", i.e. "Lines" or "Polys", which is the opendihu number of "elements", local value of rank
    global_no_t nPointsGlobal;   ///< the number of points needed for representing the mesh, global value of all rank
    global_no_t nCellsGlobal;    ///< the number of VTK "cells", i.e. "Lines" or "Polys", which is the opendihu number of "elements", global value of all ranks
    std::vector<node_no_t> nNodesLocalWithGhosts;   ///< local number of nodes including ghosts, for all dimensions

    std::vector<std::pair<std::string,int>> pointDataArrays;   ///< <name,nComponents> of PointData DataArray elements
  };
*/

  int targetDimensionality = 2;
  if (output3DMeshes)
    targetDimensionality = 3;

  VLOG(1) << "targetDimensionality: " << targetDimensionality;

  // loop over 3D or 2D meshes
  for (std::map<std::string, PolyDataPropertiesForMesh>::iterator meshPropertiesIter = meshPropertiesUnstructuredGridFile_.begin(); meshPropertiesIter != meshPropertiesUnstructuredGridFile_.end(); meshPropertiesIter++)
  {
    PolyDataPropertiesForMesh &polyDataPropertiesForMesh = meshPropertiesIter->second;

    VLOG(1) << "polyDataPropertiesForMesh A: " << polyDataPropertiesForMesh;

    if (polyDataPropertiesForMesh.dimensionality == targetDimensionality)
    {
      meshNames.insert(meshPropertiesIter->first);
      VLOG(1) << "meshName " << meshPropertiesIter->first;

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
      std::stringstream filename;
      filename << this->filenameBaseWithNo_ << "_" << targetDimensionality << "D.vtu";
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

        nCellsPreviousRanks3D_ = 0;
        nPointsPreviousRanks3D_ = 0;
        nPointsGlobal3D_ = 0;

        Control::PerformanceMeasurement::start("durationParaview3DInit");
        Control::PerformanceMeasurement::start("durationParaview3DReduction");
        MPIUtility::handleReturnValue(MPI_Exscan(&polyDataPropertiesForMesh.nCellsLocal, &nCellsPreviousRanks3D_, 1, MPI_INT, MPI_SUM, this->rankSubset_->mpiCommunicator()), "MPI_Exscan");
        MPIUtility::handleReturnValue(MPI_Exscan(&polyDataPropertiesForMesh.nPointsLocal, &nPointsPreviousRanks3D_, 1, MPI_INT, MPI_SUM, this->rankSubset_->mpiCommunicator()), "MPI_Exscan");
        MPIUtility::handleReturnValue(MPI_Reduce(&polyDataPropertiesForMesh.nPointsLocal, &nPointsGlobal3D_, 1, MPI_INT, MPI_SUM, 0, this->rankSubset_->mpiCommunicator()), "MPI_Reduce");

        Control::PerformanceMeasurement::stop("durationParaview3DReduction");
        Control::PerformanceMeasurement::stop("durationParaview3DInit");

        VLOG(1) << "nCellsLocal local: " << polyDataPropertiesForMesh.nCellsLocal << ", prefix sum: " << nCellsPreviousRanks3D_;
        VLOG(1) << "nPointsLocal local: " << polyDataPropertiesForMesh.nPointsLocal << ", prefix sum: " << nPointsPreviousRanks3D_ << ", global: " << nPointsGlobal3D_;
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
      std::vector<int> connectivityValues(nNodesPerCell*polyDataPropertiesForMesh.nCellsLocal);

      VLOG(1) << "n connectivity values from unstructured: " << polyDataPropertiesForMesh.unstructuredMeshConnectivityValues.size();
      VLOG(1) << "nCellsLocal: " << polyDataPropertiesForMesh.nCellsLocal;

      if (!polyDataPropertiesForMesh.unstructuredMeshConnectivityValues.empty())
      {
        // if connectivity values are already explicitly given, this is the case if we have an unstructured mesh to output
        assert(polyDataPropertiesForMesh.unstructuredMeshConnectivityValues.size() == connectivityValues.size());

        VLOG(1) << "connectivityValues is initialized to " << connectivityValues.size() << ", values: " << connectivityValues;
        VLOG(1) << "now copy " << polyDataPropertiesForMesh.unstructuredMeshConnectivityValues.size() << ", values: " << polyDataPropertiesForMesh.unstructuredMeshConnectivityValues;
        std::copy(polyDataPropertiesForMesh.unstructuredMeshConnectivityValues.begin(), polyDataPropertiesForMesh.unstructuredMeshConnectivityValues.end(), connectivityValues.begin());
      }
      else
      {
        // for structured meshes create connectivity values
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
                if (elementIndex*8 + 7 >= connectivityValues.size())
                {
                  LOG(FATAL) << elementIndex*8 + 7 << ">= " << connectivityValues.size() << ", connectivityValues are not large enough: " << connectivityValues.size() << ", but "
                    << nNodesLocalWithGhosts[0]-1 << "x" << nNodesLocalWithGhosts[1]-1 << "x" << nNodesLocalWithGhosts[2]-1 << " = "
                    << (nNodesLocalWithGhosts[0]-1)*(nNodesLocalWithGhosts[1]-1)*(nNodesLocalWithGhosts[2]-1) << " elements";
                }

                connectivityValues[elementIndex*8 + 0] = nPointsPreviousRanks3D_
                  + indexZ*nNodesLocalWithGhosts[0]*nNodesLocalWithGhosts[1]
                  + indexY*nNodesLocalWithGhosts[0] + indexX;
                connectivityValues[elementIndex*8 + 1] = nPointsPreviousRanks3D_
                  + indexZ*nNodesLocalWithGhosts[0]*nNodesLocalWithGhosts[1]
                  + indexY*nNodesLocalWithGhosts[0] + indexX + 1;
                connectivityValues[elementIndex*8 + 2] = nPointsPreviousRanks3D_
                  + indexZ*nNodesLocalWithGhosts[0]*nNodesLocalWithGhosts[1]
                  + (indexY+1)*nNodesLocalWithGhosts[0] + indexX + 1;
                connectivityValues[elementIndex*8 + 3] = nPointsPreviousRanks3D_
                  + indexZ*nNodesLocalWithGhosts[0]*nNodesLocalWithGhosts[1]
                  + (indexY+1)*nNodesLocalWithGhosts[0] + indexX;
                connectivityValues[elementIndex*8 + 4] = nPointsPreviousRanks3D_
                  + (indexZ+1)*nNodesLocalWithGhosts[0]*nNodesLocalWithGhosts[1]
                  + indexY*nNodesLocalWithGhosts[0] + indexX;
                connectivityValues[elementIndex*8 + 5] = nPointsPreviousRanks3D_
                  + (indexZ+1)*nNodesLocalWithGhosts[0]*nNodesLocalWithGhosts[1]
                  + indexY*nNodesLocalWithGhosts[0] + indexX + 1;
                connectivityValues[elementIndex*8 + 6] = nPointsPreviousRanks3D_
                  + (indexZ+1)*nNodesLocalWithGhosts[0]*nNodesLocalWithGhosts[1]
                  + (indexY+1)*nNodesLocalWithGhosts[0] + indexX + 1;
                connectivityValues[elementIndex*8 + 7] = nPointsPreviousRanks3D_
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
              if (elementIndex*4 + 3 >= connectivityValues.size())
              {
                LOG(FATAL) << elementIndex*4 + 3 << ">= " << connectivityValues.size() << ", connectivityValues are not large enough: " << connectivityValues.size() << ", but "
                  << nNodesLocalWithGhosts[0]-1 << "x" << nNodesLocalWithGhosts[1]-1 << " = "
                  << (nNodesLocalWithGhosts[0]-1)*(nNodesLocalWithGhosts[1]-1) << " elements";
              }

              connectivityValues[elementIndex*4 + 0] = nPointsPreviousRanks3D_
                + indexY*nNodesLocalWithGhosts[0] + indexX;
              connectivityValues[elementIndex*4 + 1] = nPointsPreviousRanks3D_
                + indexY*nNodesLocalWithGhosts[0] + indexX + 1;
              connectivityValues[elementIndex*4 + 2] = nPointsPreviousRanks3D_
                + (indexY+1)*nNodesLocalWithGhosts[0] + indexX + 1;
              connectivityValues[elementIndex*4 + 3] = nPointsPreviousRanks3D_
                + (indexY+1)*nNodesLocalWithGhosts[0] + indexX;
            }
          }
        }
      }

      VLOG(1) << "nPointsPreviousRanks3D_: " << nPointsPreviousRanks3D_ << ", nCellsPreviousRanks3D_: " << nCellsPreviousRanks3D_;
      //VLOG(1) << "nNodesLocalWithGhosts: " << nNodesLocalWithGhosts[0]-1 << "x" << nNodesLocalWithGhosts[1]-1 << "x" << nNodesLocalWithGhosts[2]-1;
      VLOG(1) << "connectivity: " << connectivityValues;

      // setup offset array
      std::vector<int> offsetValues(polyDataPropertiesForMesh.nCellsLocal);
      for (int i = 0; i < polyDataPropertiesForMesh.nCellsLocal; i++)
      {
        offsetValues[i] = nNodesPerCell*nCellsPreviousRanks3D_ + nNodesPerCell*i + nNodesPerCell;    // specifies the end, i.e. one after the last, of the last of nodes for each element
      }

      VLOG(1) << "offsetValues: " << offsetValues;

      // collect all data for the field variables, organized by field variable names
      std::map<std::string, std::vector<double>> fieldVariableValues;
      std::set<std::string> currentMeshName;
      currentMeshName.insert(meshPropertiesIter->first);
      ParaviewLoopOverTuple::loopGetNodalValues<OutputFieldVariablesType>(fieldVariables, currentMeshName, fieldVariableValues);

      // if next assertion fails, output why for debugging
      if (fieldVariableValues.size() != polyDataPropertiesForMesh.pointDataArrays.size())
      {
        LOG(DEBUG) << "n field variable values: " << fieldVariableValues.size() << ", n point data arrays: "
          << polyDataPropertiesForMesh.pointDataArrays.size();
        LOG(DEBUG) << "mesh name: " << currentMeshName;
        std::stringstream pointDataArraysNames;
        for (int i = 0; i < polyDataPropertiesForMesh.pointDataArrays.size(); i++)
        {
          pointDataArraysNames << polyDataPropertiesForMesh.pointDataArrays[i].first << " ";
        }
        LOG(DEBUG) << "pointDataArraysNames: " <<  pointDataArraysNames.str();
        LOG(DEBUG) << "OutputFieldVariablesType: " << StringUtility::demangle(typeid(OutputFieldVariablesType).name());
      }

      assert(fieldVariableValues.size() == polyDataPropertiesForMesh.pointDataArrays.size());

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
      ParaviewLoopOverTuple::loopGetGeometryFieldNodalValues<OutputFieldVariablesType>(fieldVariables, currentMeshName, geometryFieldValues);

      VLOG(1) << "currentMeshName: " << currentMeshName << ", rank " << this->rankSubset_->ownRankNo() << ", n geometryFieldValues: " << geometryFieldValues.size();
      if (geometryFieldValues.size() == 0)
      {
        LOG(FATAL) << "There is no geometry field. You have to provide a geomteryField in the field variables returned by getOutputFieldVariables!";
      }

      int nOutputFileParts = 5 + polyDataPropertiesForMesh.pointDataArrays.size();

      // create the basic structure of the output file
      std::vector<std::stringstream> outputFileParts(nOutputFileParts);
      int outputFilePartNo = 0;
      outputFileParts[outputFilePartNo] << "<?xml version=\"1.0\"?>" << std::endl
        << "<!-- " << DihuContext::versionText() << " " << DihuContext::metaText()
        << ", currentTime: " << this->currentTime_ << ", timeStepNo: " << this->timeStepNo_ << " -->" << std::endl
        << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">" << std::endl    // intel cpus are LittleEndian
        << std::string(1, '\t') << "<UnstructuredGrid>" << std::endl;

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
    else
    {
      VLOG(1) << "skip mesh " << meshPropertiesIter->first << " because " << polyDataPropertiesForMesh.dimensionality << " != " << targetDimensionality;

    }
  }
}

} // namespace
