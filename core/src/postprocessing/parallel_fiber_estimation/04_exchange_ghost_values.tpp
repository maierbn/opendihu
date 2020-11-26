#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
bool ParallelFiberEstimation<BasisFunctionType>::
neighbourExists(const std::array<bool,4> &subdomainIsAtBorder, Mesh::face_or_edge_t faceOrEdge)
{
  // determine if there is a neighboring rank

  // if it is a face
  if (faceOrEdge < 4)
  {
    if (!subdomainIsAtBorder[faceOrEdge])
    {
      return true;
    }
  }
  else
  {
    // it is an edge
    if (
         (faceOrEdge == Mesh::face_or_edge_t::edge0Minus1Minus && !subdomainIsAtBorder[Mesh::face_t::face0Minus] && !subdomainIsAtBorder[Mesh::face_t::face1Minus])
      || (faceOrEdge == Mesh::face_or_edge_t::edge0Plus1Minus  && !subdomainIsAtBorder[Mesh::face_t::face0Plus]  && !subdomainIsAtBorder[Mesh::face_t::face1Minus])
      || (faceOrEdge == Mesh::face_or_edge_t::edge0Minus1Plus  && !subdomainIsAtBorder[Mesh::face_t::face0Minus] && !subdomainIsAtBorder[Mesh::face_t::face1Plus])
      || (faceOrEdge == Mesh::face_or_edge_t::edge0Plus1Plus   && !subdomainIsAtBorder[Mesh::face_t::face0Plus]  && !subdomainIsAtBorder[Mesh::face_t::face1Plus])
    )
    {
      return true;
    }
  }
  return false;
}

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
printRanksInNeighbourhood(const std::array<bool,4> &subdomainIsAtBorder)
{
  std::stringstream s;

  // output neighbouring ranks and own rank for debugging
  // top row
  if (subdomainIsAtBorder[Mesh::face_t::face1Plus])
  {
    // left
    if (subdomainIsAtBorder[Mesh::face_t::face0Minus])
    {
      s << "+--";
    }
    else
    {
      s << "---";
    }

    // center
    s << "---";

    // right
    if (subdomainIsAtBorder[Mesh::face_t::face0Plus])
    {
      s << "+";
    }
    else
    {
      s << "---";
    }
  }
  else
  {
    // left
    if (subdomainIsAtBorder[Mesh::face_t::face0Minus])
    {
      s << "|  ";
    }
    else
    {
      int neighbourRankNo = meshPartition_->neighbourRank(Mesh::face_or_edge_t::edge0Minus1Plus);
      s << " " << std::setw(2) << neighbourRankNo;
    }

    // center
    int neighbourRankNo = meshPartition_->neighbourRank(Mesh::face_or_edge_t::faceEdge1Plus);
    s << " " << std::setw(2) << neighbourRankNo;

    // right
    if (subdomainIsAtBorder[Mesh::face_t::face0Plus])
    {
      s << "|";
    }
    else
    {
      int neighbourRankNo = meshPartition_->neighbourRank(Mesh::face_or_edge_t::edge0Plus1Plus);
      s << " " << std::setw(2) << neighbourRankNo;
    }
  }

  s << "\n";

  // middle row
  // left
  if (subdomainIsAtBorder[Mesh::face_t::face0Minus])
  {
    s << "|  ";
  }
  else
  {
    int neighbourRankNo = meshPartition_->neighbourRank(Mesh::face_or_edge_t::faceEdge0Minus);
    s << " " << std::setw(2) << neighbourRankNo;
  }

  // center
  s << " " << std::setw(2) << currentRankSubset_->ownRankNo();

  // right
  if (subdomainIsAtBorder[Mesh::face_t::face0Plus])
  {
    s << "|";
  }
  else
  {
    int neighbourRankNo = meshPartition_->neighbourRank(Mesh::face_or_edge_t::faceEdge0Plus);
    s << " " << std::setw(2) << neighbourRankNo;
  }

  s << "\n";

  // bottom row
  if (subdomainIsAtBorder[Mesh::face_t::face1Minus])
  {
    // left
    if (subdomainIsAtBorder[Mesh::face_t::face0Minus])
    {
      s << "+--";
    }
    else
    {
      s << "---";
    }

    // center
    s << "---";

    // right
    if (subdomainIsAtBorder[Mesh::face_t::face0Plus])
    {
      s << "+";
    }
    else
    {
      s << "---";
    }
  }
  else
  {
    // left
    if (subdomainIsAtBorder[Mesh::face_t::face0Minus])
    {
      s << "|  ";
    }
    else
    {
      int neighbourRankNo = meshPartition_->neighbourRank(Mesh::face_or_edge_t::edge0Minus1Minus);
      s << " " << std::setw(2) << neighbourRankNo;
    }

    // center
    int neighbourRankNo = meshPartition_->neighbourRank(Mesh::face_or_edge_t::faceEdge1Minus);
    s << " " << std::setw(2) << neighbourRankNo;

    // right
    if (subdomainIsAtBorder[Mesh::face_t::face0Plus])
    {
      s << "|";
    }
    else
    {
      int neighbourRankNo = meshPartition_->neighbourRank(Mesh::face_or_edge_t::edge0Plus1Minus);
      s << " " << std::setw(2) << neighbourRankNo;
    }
  }
  LOG(DEBUG) << "current ranks in neighbourhood:\n" << s.str();
}

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
exchangeGhostValues(const std::array<bool,4> &subdomainIsAtBorder)
{
  struct GhostValues
  {
    std::vector<double> nodePositionValues;  //< the values of the node positions of the ghost elements, in sequential order like [x,x,x,x, ... y,y,y,y ... z,z,z,...]
    std::vector<double> solutionValues;      //< values of solution field variable inside the ghost elements
    std::vector<double> gradientValues;      //< values of the gradient field variable, consecutive like nodePositionValues
    std::array<element_no_t,3> nElementsPerCoordinateDirection;   //< size of the ghost mesh
  };
  std::array<GhostValues,8> ghostValuesBuffer;  //< [faceIndex], data for meshes containing ghost elements for the sides, face0Minus, face0Plus, face1Minus, face1Plus
  std::array<GhostValues,8> boundaryValues;     //< [faceIndex], data to be send

  std::vector<MPI_Request> receiveRequests;
  std::vector<MPI_Request> sendRequests;

  std::vector<Mesh::face_or_edge_t> faces = {
    Mesh::face_or_edge_t::faceEdge0Minus, Mesh::face_or_edge_t::faceEdge0Plus, 
    Mesh::face_or_edge_t::faceEdge1Minus, Mesh::face_or_edge_t::faceEdge1Plus,
    Mesh::face_or_edge_t::edge0Minus1Minus, Mesh::face_or_edge_t::edge0Plus1Plus,
    Mesh::face_or_edge_t::edge0Plus1Minus, Mesh::face_or_edge_t::edge0Minus1Plus
  };

  // communicate ghost elements to neighbouring subdomains
  // determine elements on own domain that are to be send to the neighbouring domain
  // loop over faces: face0Minus = 0, face0Plus, face1Minus, face1Plus
  // and edges: edge0Minus1Minus, edge0Plus1Minus, edge0Minus1Plus, edge0Plus1Plus
  for (int faceIndex = 0; faceIndex != 8; faceIndex++)
  {
    Mesh::face_or_edge_t faceOrEdge = faces[faceIndex];

    if (neighbourExists(subdomainIsAtBorder, faceOrEdge))
    {
      LOG(DEBUG) << "determine ghost elements for face " << Mesh::getString((Mesh::face_or_edge_t)faceOrEdge);

      // get information about neighbouring rank and boundary elements for face
      int neighbourRankNo;
      std::vector<dof_no_t> dofNos;
      meshPartition_->getBoundaryElements((Mesh::face_or_edge_t)faceOrEdge, ghostLayerWidth_, neighbourRankNo, ghostValuesBuffer[faceIndex].nElementsPerCoordinateDirection, dofNos);

      LOG(DEBUG) << "getBoundaryElements for face " << Mesh::getString((Mesh::face_or_edge_t)faceOrEdge) << " returned neighbourRankNo=" << neighbourRankNo
        << ", nElementsPerCoordinateDirection: " << ghostValuesBuffer[faceIndex].nElementsPerCoordinateDirection << ", " << dofNos.size() << " dofs: "
        << dofNos;

      // get relevant values in own domain that will be send to the neighbouring domain
      problem_->data().functionSpace()->geometryField().getValues(dofNos, boundaryValues[faceIndex].nodePositionValues);
      problem_->data().solution()->getValues(dofNos, boundaryValues[faceIndex].solutionValues);
      data_.gradient()->getValues(dofNos, boundaryValues[faceIndex].gradientValues);
    }
  }

  // this barrier is needed such that VecGhostUpdate does not interfere with the following asynchronous communication
  MPI_Barrier(this->currentRankSubset_->mpiCommunicator());
  LOG(DEBUG) << "determined boundary elements, now communicate";

  // for debugging, output neighboring ranks
  printRanksInNeighbourhood(subdomainIsAtBorder);

  // blocking communication
#if 0
  // left: recv-send, right: send-recv
  // bottom: recv-send, top: send-recv
  // bottom-left: recv-send, top-right: send-recv
  // bottom-right: recv-send, top-left: send-recv
  // loop over the faces and edges
  for (int faceIndex = 0; faceIndex != 8; faceIndex++)
  {
    Mesh::face_or_edge_t faceOrEdge = faces[faceIndex];

    if (neighbourExists(subdomainIsAtBorder, faceOrEdge))
    {
      int neighbourRankNo = meshPartition_->neighbourRank((Mesh::face_or_edge_t)faceOrEdge);
      assert (neighbourRankNo != -1);

      int nNodePositionValues = boundaryValues[faceIndex].nodePositionValues.size();
      int nSolutionValues = boundaryValues[faceIndex].solutionValues.size();
      int nGradientValues = boundaryValues[faceIndex].gradientValues.size();
      assert(nSolutionValues*3 == nGradientValues);
      assert(nNodePositionValues == nGradientValues);
      int tag = 0;

      // checkerboard pattern for the own domain
      for (int i = 0; i != 2; i++)
      {
        if (faceIndex % 2 == i)
        {
          // receive from neighbouring process
          // blocking receive call to receive node position values
          ghostValuesBuffer[faceIndex].nodePositionValues.resize(nNodePositionValues);
          tag = currentRankSubset_->ownRankNo()*100 + neighbourRankNo*10000 + level_*10 + 1;
          LOG(DEBUG) << "receive " << nNodePositionValues << " (" << ghostValuesBuffer[faceIndex].nodePositionValues.size() << ") from rank " << neighbourRankNo << " (" << Mesh::getString(faceOrEdge) << ") (tag: " << tag << ")";
          MPIUtility::handleReturnValue(MPI_Recv(ghostValuesBuffer[faceIndex].nodePositionValues.data(), nNodePositionValues, MPI_DOUBLE,
                                                  neighbourRankNo, tag, currentRankSubset_->mpiCommunicator(), MPI_STATUS_IGNORE), "MPI_Recv");


          // blocking receive call to receive solution values
          ghostValuesBuffer[faceIndex].solutionValues.resize(nSolutionValues);
          tag = currentRankSubset_->ownRankNo()*100 + neighbourRankNo*10000 + level_*10 + 2;
          LOG(DEBUG) << "receive " << nSolutionValues << " (" << ghostValuesBuffer[faceIndex].solutionValues.size() << ") from rank " << neighbourRankNo << " (" << Mesh::getString(faceOrEdge) << ") (tag: " << tag << ")";
          MPIUtility::handleReturnValue(MPI_Recv(ghostValuesBuffer[faceIndex].solutionValues.data(), nSolutionValues, MPI_DOUBLE,
                                                  neighbourRankNo, tag, currentRankSubset_->mpiCommunicator(), MPI_STATUS_IGNORE), "MPI_Recv");


          // blocking receive call to receive gradient values
          ghostValuesBuffer[faceIndex].gradientValues.resize(nGradientValues);
          tag = currentRankSubset_->ownRankNo()*100 + neighbourRankNo*10000 + level_*10 + 3;
          LOG(DEBUG) << "receive " << nGradientValues << " (" << ghostValuesBuffer[faceIndex].gradientValues.size() << ") from rank " << neighbourRankNo << " (" << Mesh::getString(faceOrEdge) << ") (tag: " << tag << ")";
          MPIUtility::handleReturnValue(MPI_Recv(ghostValuesBuffer[faceIndex].gradientValues.data(), nGradientValues, MPI_DOUBLE,
                                                  neighbourRankNo, tag, currentRankSubset_->mpiCommunicator(), MPI_STATUS_IGNORE), "MPI_Recv");

          LOG(DEBUG) << "receive from rank " << neighbourRankNo << " completed";
        }
        else
        {
          // send values to neighbouring process
          // blocking send call to send solution values
          tag = currentRankSubset_->ownRankNo()*10000 + neighbourRankNo*100 + level_*10 + 1;
          LOG(DEBUG) << "send " << nNodePositionValues << " (" << boundaryValues[faceIndex].nodePositionValues.size() << ") to rank " << neighbourRankNo << " (" << Mesh::getString(faceOrEdge) << ") (tag: " << tag << ")";
          MPIUtility::handleReturnValue(MPI_Send(boundaryValues[faceIndex].nodePositionValues.data(), nNodePositionValues, MPI_DOUBLE,
                                                  neighbourRankNo, tag, currentRankSubset_->mpiCommunicator()), "MPI_Send");


          // blocking send call to send solution values
          tag = currentRankSubset_->ownRankNo()*10000 + neighbourRankNo*100 + level_*10 + 2;
          LOG(DEBUG) << "send " << nSolutionValues << " (" << boundaryValues[faceIndex].solutionValues.size() << ") to rank " << neighbourRankNo << " (" << Mesh::getString(faceOrEdge) << ") (tag: " << tag << ")";
          MPIUtility::handleReturnValue(MPI_Send(boundaryValues[faceIndex].solutionValues.data(), nSolutionValues, MPI_DOUBLE,
                                                  neighbourRankNo, tag, currentRankSubset_->mpiCommunicator()), "MPI_Send");


          // blocking send call to send gradient values
          tag = currentRankSubset_->ownRankNo()*10000 + neighbourRankNo*100 + level_*10 + 3;
          LOG(DEBUG) << "send " << nGradientValues << " (" << boundaryValues[faceIndex].gradientValues.size() << ") to rank " << neighbourRankNo << " (" << Mesh::getString(faceOrEdge) << ") (tag: " << tag << ")";
          MPIUtility::handleReturnValue(MPI_Send(boundaryValues[faceIndex].gradientValues.data(), nGradientValues, MPI_DOUBLE,
                                                  neighbourRankNo, tag, currentRankSubset_->mpiCommunicator()), "MPI_Send");

          LOG(DEBUG) << "send to rank " << neighbourRankNo << " completed";
        }
      }
    }
  }
#endif 

// non-blocking communication, does not work with 64 processes
#if 1
  sendRequests.clear();
  receiveRequests.clear();

  // loop over faces and communicate ghost elements to the neighbouring ranks
  for (int faceIndex = 0; faceIndex != 8; faceIndex++)
  {
    Mesh::face_or_edge_t faceOrEdge = faces[faceIndex];

    if (neighbourExists(subdomainIsAtBorder, faceOrEdge))
    {
      int neighbourRankNo = meshPartition_->neighbourRank((Mesh::face_or_edge_t)faceOrEdge);
      assert (neighbourRankNo != -1);

#if 0
      // output sent ghost elements for debugging
      std::stringstream s;
      s << "04_ghost_elements_face_" << Mesh::getString((Mesh::face_or_edge_t)faceOrEdge);
      PyObject_CallFunction(functionOutputGhostElements_, "s i i O O f", s.str().c_str(), currentRankSubset_->ownRankNo(), level_,
                            PythonUtility::convertToPython<std::vector<double>>::get(boundaryValues[faceIndex].nodePositionValues),
                            PythonUtility::convertToPython<std::array<element_no_t,3>>::get(ghostValuesBuffer[faceIndex].nElementsPerCoordinateDirection), 0.05);
      PythonUtility::checkForError();

#endif

      int nNodePositionValues = boundaryValues[faceIndex].nodePositionValues.size();
      int nSolutionValues = boundaryValues[faceIndex].solutionValues.size();
      int nGradientValues = boundaryValues[faceIndex].gradientValues.size();
      assert(nSolutionValues*3 == nGradientValues);
      assert(nNodePositionValues == nGradientValues);

      LOG(DEBUG) << "exchange ghosts with neighbour " << neighbourRankNo << " (" << Mesh::getString(faceOrEdge) << ") : nNodePositionValues=" << nNodePositionValues << ", nSolutionValues=" << nSolutionValues << ", nGradientValues=" << nGradientValues;

      // if no checkpoint is used here for debugging
#if !defined(USE_CHECKPOINT_GHOST_MESH)

      LOG(DEBUG) << "receive from rank " << neighbourRankNo;

      // receive from neighbouring process
      // post non-blocking receive call to receive node position values
      MPI_Request receiveRequest;
      ghostValuesBuffer[faceIndex].nodePositionValues.resize(nNodePositionValues);
      MPIUtility::handleReturnValue(MPI_Irecv(ghostValuesBuffer[faceIndex].nodePositionValues.data(), nNodePositionValues, MPI_DOUBLE,
                                              neighbourRankNo, 0, currentRankSubset_->mpiCommunicator(), &receiveRequest), "MPI_Irecv");
      receiveRequests.push_back(receiveRequest);

      // post non-blocking receive call to receive solution values
      ghostValuesBuffer[faceIndex].solutionValues.resize(nSolutionValues);
      MPIUtility::handleReturnValue(MPI_Irecv(ghostValuesBuffer[faceIndex].solutionValues.data(), nSolutionValues, MPI_DOUBLE,
                                              neighbourRankNo, 0, currentRankSubset_->mpiCommunicator(), &receiveRequest), "MPI_Irecv");
      receiveRequests.push_back(receiveRequest);

      // post non-blocking receive call to receive gradient values
      ghostValuesBuffer[faceIndex].gradientValues.resize(nGradientValues);
      MPIUtility::handleReturnValue(MPI_Irecv(ghostValuesBuffer[faceIndex].gradientValues.data(), nGradientValues, MPI_DOUBLE,
                                              neighbourRankNo, 0, currentRankSubset_->mpiCommunicator(), &receiveRequest), "MPI_Irecv");
      receiveRequests.push_back(receiveRequest);

      LOG(DEBUG) << "receive from rank " << neighbourRankNo << " (" << Mesh::getString(faceOrEdge) << ") completed";

      LOG(DEBUG) << "send to rank " << neighbourRankNo << " (" << Mesh::getString(faceOrEdge) << ")";

      // send values to neighbouring process
      // post non-blocking send call to send solution values
      MPI_Request sendRequest;
      MPIUtility::handleReturnValue(MPI_Isend(boundaryValues[faceIndex].nodePositionValues.data(), nNodePositionValues, MPI_DOUBLE,
                                              neighbourRankNo, 0, currentRankSubset_->mpiCommunicator(), &sendRequest), "MPI_Isend");
      sendRequests.push_back(sendRequest);

      // post non-blocking send call to send solution values
      MPIUtility::handleReturnValue(MPI_Isend(boundaryValues[faceIndex].solutionValues.data(), nSolutionValues, MPI_DOUBLE,
                                              neighbourRankNo, 0, currentRankSubset_->mpiCommunicator(), &sendRequest), "MPI_Isend");
      sendRequests.push_back(sendRequest);

      // post non-blocking send call to send gradient values
      MPIUtility::handleReturnValue(MPI_Isend(boundaryValues[faceIndex].gradientValues.data(), nGradientValues, MPI_DOUBLE,
                                              neighbourRankNo, 0, currentRankSubset_->mpiCommunicator(), &sendRequest), "MPI_Isend");
      sendRequests.push_back(sendRequest);

      LOG(DEBUG) << "send to rank " << neighbourRankNo << " (" << Mesh::getString(faceOrEdge) << ") completed";
#endif
    }
  }

  // wait for non-blocking communication to finish
  MPIUtility::handleReturnValue(MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");
  MPIUtility::handleReturnValue(MPI_Waitall(receiveRequests.size(), receiveRequests.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");

  LOG(DEBUG) << "waitall (" << sendRequests.size() << " send requests, " << receiveRequests.size() << " receiveRefquests) complete";
#endif
  LOG(DEBUG) << "ghost exchange communication done, currentRankSubset_: " << *currentRankSubset_;


  //MPI_Barrier(currentRankSubset_->mpiCommunicator());
  //LOG(FATAL) << "end after get ghost elements";

  // for debugging, write and load checkpoints and output data for debugging

  // loop over the faces and edges
  for (int faceIndex = 0; faceIndex != 8; faceIndex++)
  {
    Mesh::face_or_edge_t faceOrEdge = faces[faceIndex];

    if (neighbourExists(subdomainIsAtBorder, faceOrEdge))
    {
#if defined(WRITE_CHECKPOINT_GHOST_MESH) || defined(USE_CHECKPOINT_GHOST_MESH) || !defined(NDEBUG)
      int neighbourRankNo = meshPartition_->neighbourRank((Mesh::face_or_edge_t)faceOrEdge);
#endif

#ifdef WRITE_CHECKPOINT_GHOST_MESH
      // write checkpoint of all data at this point in the code (for debugging)
      std::stringstream filenameOut;
      filenameOut << "checkpoints/checkpoint_ghost_mesh_l" << level_ << "_" << neighbourRankNo << "_" << Mesh::getString((Mesh::face_or_edge_t)faceOrEdge) << ".txt";
      std::ofstream fileOut;
      fileOut.open(filenameOut.str().c_str(), std::ios::out | std::ios::trunc);

      assert(fileOut.is_open());

      fileOut << ghostValuesBuffer[faceIndex].nodePositionValues.size() << " "
        << ghostValuesBuffer[faceIndex].solutionValues.size() << " "
        << ghostValuesBuffer[faceIndex].gradientValues.size() << " ";

      LOG(DEBUG) << " output " << ghostValuesBuffer[faceIndex].nodePositionValues.size() << " node position values, "
        << ghostValuesBuffer[faceIndex].solutionValues.size() << " solution values, " << ghostValuesBuffer[faceIndex].gradientValues.size()
        << " gradient values to file " << filenameOut.str();
      LOG(DEBUG) << "ghostValuesBuffer[faceIndex].nodePositionValues: " << ghostValuesBuffer[faceIndex].nodePositionValues;
      LOG(DEBUG) << "ghostValuesBuffer[faceIndex].solutionValues: " << ghostValuesBuffer[faceIndex].solutionValues;
      LOG(DEBUG) << "ghostValuesBuffer[faceIndex].gradientValues: " << ghostValuesBuffer[faceIndex].gradientValues;

      for (int i = 0; i < ghostValuesBuffer[faceIndex].nodePositionValues.size(); i++)
      {
        if (i % int(ghostValuesBuffer[faceIndex].nodePositionValues.size()/3) == 0)
          fileOut << std::endl;

        fileOut << ghostValuesBuffer[faceIndex].nodePositionValues[i] << " ";
      }
      fileOut << std::endl;

      for (int i = 0; i < ghostValuesBuffer[faceIndex].solutionValues.size(); i++)
        fileOut << ghostValuesBuffer[faceIndex].solutionValues[i] << " ";

      fileOut << std::endl;

      for (int i = 0; i < ghostValuesBuffer[faceIndex].gradientValues.size(); i++)
        fileOut << ghostValuesBuffer[faceIndex].gradientValues[i] << " ";

      fileOut << std::endl;

      fileOut.close();
      LOG(DEBUG) << "Wrote ghost meshes (" << ghostValuesBuffer[faceIndex].nodePositionValues.size() << " node position values, "
        << ghostValuesBuffer[faceIndex].solutionValues.size() << " solution values, " << ghostValuesBuffer[faceIndex].gradientValues.size()
        << " gradient values) to checkpoint file " << filenameOut.str();
#endif

#ifdef USE_CHECKPOINT_GHOST_MESH
      // load all data from a checkpoint as this point in the code (for debugging)
      std::stringstream filenameIn;
      filenameIn << "checkpoints/checkpoint_ghost_mesh_l" << level_ << "_" << neighbourRankNo << "_" << Mesh::getString((Mesh::face_or_edge_t)faceOrEdge) << ".txt";
      std::ifstream fileIn;
      fileIn.open(filenameIn.str().c_str(), std::ios::in);

      assert(fileIn.is_open());

      int size1, size2, size3;
      fileIn >> size1 >> size2 >> size3;
      ghostValuesBuffer[faceIndex].nodePositionValues.resize(size1);
      ghostValuesBuffer[faceIndex].solutionValues.resize(size2);
      ghostValuesBuffer[faceIndex].gradientValues.resize(size3);

      for (int i = 0; i < ghostValuesBuffer[faceIndex].nodePositionValues.size(); i++)
        fileIn >> ghostValuesBuffer[faceIndex].nodePositionValues[i];

      for (int i = 0; i < ghostValuesBuffer[faceIndex].solutionValues.size(); i++)
        fileIn >> ghostValuesBuffer[faceIndex].solutionValues[i];

      for (int i = 0; i < ghostValuesBuffer[faceIndex].gradientValues.size(); i++)
        fileIn >> ghostValuesBuffer[faceIndex].gradientValues[i];

      fileIn.close();
      LOG(DEBUG) << "Loaded ghost meshes (" << ghostValuesBuffer[faceIndex].nodePositionValues.size() << " node position values, "
        << ghostValuesBuffer[faceIndex].solutionValues.size() << " solution values, " << ghostValuesBuffer[faceIndex].gradientValues.size()
        << " gradient values) from checkpoint file " << filenameIn.str();


      int elementNo = 392;
      // element coordinates 0,0,49
      if (neighbourRankNo == 1 && Mesh::getString((Mesh::face_or_edge_t)faceOrEdge) == "0+")
      {
        LOG(DEBUG) << "special element " << elementNo;
        //int
        // 0+
        // element coordinates: (0,0,49) / (1,8,50)
        int dof[8] = {49*(2*9), 49*(2*9)+1, 49*(2*9) + 1*2, 49*(2*9) + 1*2 + 1,
                      50*(2*9), 50*(2*9)+1, 50*(2*9) + 1*2, 50*(2*9) + 1*2 + 1};
        for (int i = 0; i < 8; i++)
        {
          LOG(DEBUG) << ghostValuesBuffer[faceIndex].nodePositionValues[dof[i]] << ","
            << ghostValuesBuffer[faceIndex].nodePositionValues[918 + dof[i]] << ","
            << ghostValuesBuffer[faceIndex].nodePositionValues[918*2 + dof[i]];
        }
      }

#endif

      // for debug target
#ifndef NDEBUG
      // output received ghost elements
      if (neighbourRankNo != -1)
      {
        std::stringstream s;
        s << "05_received_ghost_elements_face_" << Mesh::getString((Mesh::face_or_edge_t)faceOrEdge);
        PyObject_CallFunction(functionOutputGhostElements_, "s i i O O f", s.str().c_str(), currentRankSubset_->ownRankNo(), level_,
                              PythonUtility::convertToPython<std::vector<double>>::get(ghostValuesBuffer[faceIndex].nodePositionValues),
                              PythonUtility::convertToPython<std::array<element_no_t,3>>::get(ghostValuesBuffer[faceIndex].nElementsPerCoordinateDirection), 0.05);
        PythonUtility::checkForError();
      }
#endif
    }
  }

  // handle received values
  // create ghost element meshes from received data

  // create a new communicator that contains only the current rank, because the ghost mesh is not partitioned
  LOG(DEBUG) << "split rank subset, own rank: " << currentRankSubset_->ownRankNo();
  std::shared_ptr<Partition::RankSubset> rankSubsetSingleRank = std::make_shared<Partition::RankSubset>(currentRankSubset_->ownRankNo(), currentRankSubset_);

  // loop over the faces and edges
  for (int faceIndex = 0; faceIndex != 8; faceIndex++)
  {
    Mesh::face_or_edge_t faceOrEdge = faces[faceIndex];

    std::shared_ptr<FunctionSpaceType> ghostMesh;

    // if values were received
    if (!ghostValuesBuffer[faceIndex].nodePositionValues.empty())
    {
      // compose name for new ghost mesh
      std::stringstream meshName;
      meshName << this->functionSpace_->meshName() << "_ghost_" << Mesh::getString((Mesh::face_or_edge_t)faceOrEdge);

      // transform the node position data from vector of double to vector of Vec3
      int nNodes = ghostValuesBuffer[faceIndex].nodePositionValues.size() / 3;
      std::vector<Vec3> nodePositions(nNodes);

      for (node_no_t nodeIndex = 0; nodeIndex < nNodes; nodeIndex++)
      {
        nodePositions[nodeIndex][0] = ghostValuesBuffer[faceIndex].nodePositionValues[0*nNodes + nodeIndex];
        nodePositions[nodeIndex][1] = ghostValuesBuffer[faceIndex].nodePositionValues[1*nNodes + nodeIndex];
        nodePositions[nodeIndex][2] = ghostValuesBuffer[faceIndex].nodePositionValues[2*nNodes + nodeIndex];
      }

      // the number of elements per coordinate direction was set earlier by getBoundaryElements
      LOG(DEBUG) << "create ghost mesh with nElementsPerCoordinateDirection: " << ghostValuesBuffer[faceIndex].nElementsPerCoordinateDirection;

      // create ghost mesh
      std::array<int,3> nRanks({1,1,1});
      context_.partitionManager()->setRankSubsetForNextCreatedPartitioning(rankSubsetSingleRank);
      ghostMesh = context_.meshManager()->template createFunctionSpace<FunctionSpaceType>(
        meshName.str(), nodePositions, ghostValuesBuffer[faceIndex].nElementsPerCoordinateDirection, nRanks);

      LOG(DEBUG) << "created ghost mesh for face " << Mesh::getString((Mesh::face_or_edge_t)faceOrEdge) << ".";

      // create solution field variable for ghost mesh
      LOG(DEBUG) << "create solution field variable on ghost mesh";
      this->ghostMeshSolution_[faceOrEdge] = ghostMesh->template createFieldVariable<1>("solution");
      this->ghostMeshSolution_[faceOrEdge]->setValuesWithGhosts(ghostValuesBuffer[faceIndex].solutionValues);

      // create gradient field variable for ghost mesh
      LOG(DEBUG) << "create gradient field variable on ghost mesh";
      this->ghostMeshGradient_[faceOrEdge] = ghostMesh->template createFieldVariable<3>("gradient");
      int nGradientValues = ghostValuesBuffer[faceIndex].gradientValues.size()/3;

      // convert gradient values from vector of doubles to vector of Vec3's
      std::vector<Vec3> gradientValues(nGradientValues);
      for (int i = 0; i < nGradientValues; i++)
      {
        gradientValues[i][0] = ghostValuesBuffer[faceIndex].gradientValues[3*i + 0];
        gradientValues[i][1] = ghostValuesBuffer[faceIndex].gradientValues[3*i + 1];
        gradientValues[i][2] = ghostValuesBuffer[faceIndex].gradientValues[3*i + 2];
      }
      this->ghostMeshGradient_[faceOrEdge]->setValuesWithGhosts(gradientValues);
    }
    else
    {
      ghostMesh = nullptr;
      LOG(DEBUG) << "ghost mesh for face " << Mesh::getString((Mesh::face_or_edge_t)faceOrEdge) << " is null, because ghostValuesBuffer[faceIndex].nodePositionValues is empty()";
    }

    this->functionSpace_->setGhostMesh((Mesh::face_or_edge_t)faceOrEdge, ghostMesh);
    //LOG(DEBUG) << "after settings ghost mesh for face " << Mesh::getString((Mesh::face_or_edge_t)faceOrEdge) << ": ";
    //this->functionSpace_->debugOutputGhostMeshSet();
  }

}



} // namespace
