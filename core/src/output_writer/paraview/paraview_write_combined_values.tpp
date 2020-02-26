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

} // namespace
