#include "specialized_solver/fast_monodomain_solver/fast_monodomain_solver_base.h"

#include "partition/rank_subset.h"
#include "control/stimulation_logging.h"

template<int nStates, int nIntermediates>
FastMonodomainSolverBase<nStates,nIntermediates>::
FastMonodomainSolverBase(const DihuContext &context) :
  specificSettings_(context.getPythonConfig()), nestedSolvers_(context)
{
  // initialize output writers
  this->outputWriterManager_.initialize(context, specificSettings_);
}

template<int nStates, int nIntermediates>
void FastMonodomainSolverBase<nStates,nIntermediates>::
initialize()
{
  nestedSolvers_.initialize();

  // initialize motor unit numbers and firing times
  fiberDistributionFilename_ = specificSettings_.getOptionString("fiberDistributionFile", "");
  firingTimesFilename_ = specificSettings_.getOptionString("firingTimesFile", "");


  // output warning if there are output writers
  if (this->outputWriterManager_.hasOutputWriters())
  {
    LOG(WARNING) << "The FastMonodomainSolver has output writers. Those will output all states of the 0D problem every splitting step. "
      << "This may be slow! Maybe you do not want to have this because the 1D fiber output writers can also be called.";
  }

  LOG(DEBUG) << "config: " << specificSettings_;
  LOG(DEBUG) << "fiberDistributionFilename: " << fiberDistributionFilename_;
  LOG(DEBUG) << "firingTimesFilename: " << firingTimesFilename_;

  // parse firingTimesFilename_
  std::shared_ptr<Partition::RankSubset> rankSubset = nestedSolvers_.data().functionSpace()->meshPartition()->rankSubset();
  std::string firingTimesFileContents = MPIUtility::loadFile(firingTimesFilename_, rankSubset->mpiCommunicator());

  // parse file contents
  while (!firingTimesFileContents.empty())
  {
    // extract line
    std::size_t lineEndPos = firingTimesFileContents.find("\n");
    if (lineEndPos == std::string::npos)
      break;

    std::string line = firingTimesFileContents.substr(0, lineEndPos);
    firingTimesFileContents.erase(0, lineEndPos+1);

    firingEvents_.push_back(std::vector<bool>());

    // parse line
    while (!line.empty())
    {
      int entry = atoi(line.c_str());

      firingEvents_.back().push_back((bool)(entry));

      std::size_t pos = line.find_first_of("\t ");
      if (pos == std::string::npos)
        break;
      line.erase(0, pos+1);

      while (!isdigit(line[0]) && !line.empty())
        line.erase(0,1);
    }
  }

  // parse fiberDistributionFile
  std::string fiberDistributionFileContents = MPIUtility::loadFile(fiberDistributionFilename_, rankSubset->mpiCommunicator());

  // parse file contents
  while (!fiberDistributionFileContents.empty())
  {
    int motorUnitNo = atoi(fiberDistributionFileContents.c_str());
    motorUnitNo_.push_back(motorUnitNo);

    std::size_t pos = fiberDistributionFileContents.find_first_of("\t ");

    if (pos == std::string::npos)
      break;

    fiberDistributionFileContents.erase(0, pos+1);
  }

  LOG(DEBUG) << "firingEvents.size: " << firingEvents_.size();
  LOG(DEBUG) << "firingEvents_:" << firingEvents_;
  LOG(DEBUG) << "motorUnitNo_: " << motorUnitNo_;

  if (motorUnitNo_.empty())
    LOG(FATAL) << "Could not parse motor units.";

  if (firingEvents_.empty())
    LOG(FATAL) << "Could not parse firing times.";

  // initialize data structures
  std::vector<typename NestedSolversType::TimeSteppingSchemeType> &instances = nestedSolvers_.instancesLocal();

  // determine number of fibers to compute on the current rank
  int nFibers = 0;
  int fiberNo = 0;
  nFibersToCompute_ = 0;
  for (int i = 0; i < instances.size(); i++)
  {
    std::vector<TimeSteppingScheme::Heun<CellmlAdapterType>> &innerInstances
      = instances[i].timeStepping1().instancesLocal();  // TimeSteppingScheme::Heun<CellmlAdapter...
    nFibers += innerInstances.size();

    for (int j = 0; j < innerInstances.size(); j++, fiberNo++)
    {
      std::shared_ptr<FiberFunctionSpace> fiberFunctionSpace = innerInstances[j].data().functionSpace();
      std::shared_ptr<Partition::RankSubset> rankSubset = fiberFunctionSpace->meshPartition()->rankSubset();
      int computingRank = fiberNo % rankSubset->size();

      if (computingRank == rankSubset->ownRankNo())
      {
        LOG(DEBUG) << "compute (i,j)=(" << i << "," << j << "), computingRank " << computingRank 
          << ", fiberNo: " << fiberNo << ", rankSubset->size(): " << rankSubset->size() << ", own: " << rankSubset->ownRankNo();
        nFibersToCompute_++;
      }
    }
  }

  fiberData_.resize(nFibersToCompute_);
  LOG(DEBUG) << "nFibers: " << nFibers << ", nFibersToCompute_: " << nFibersToCompute_;

  // determine total number of Hodgkin-Huxley instances to compute on this rank
  nInstancesToCompute_ = 0;
  int fiberDataNo = 0;
  fiberNo = 0;
  for (int i = 0; i < instances.size(); i++)
  {
    std::vector<TimeSteppingScheme::Heun<CellmlAdapterType>> &innerInstances
      = instances[i].timeStepping1().instancesLocal();  // TimeSteppingScheme::Heun<CellmlAdapter...

    for (int j = 0; j < innerInstances.size(); j++, fiberNo++)
    {
      CellmlAdapterType &cellmlAdapter = innerInstances[j].discretizableInTime();
      int fiberNoGlobal = PythonUtility::convertFromPython<int>::get(cellmlAdapter.pySetFunctionAdditionalParameter_);

      // initialize filename in stimulation logging class from current settings
      Control::StimulationLogging logging(cellmlAdapter.specificSettings_);

      std::shared_ptr<FiberFunctionSpace> fiberFunctionSpace = innerInstances[j].data().functionSpace();

      std::shared_ptr<Partition::RankSubset> rankSubset = fiberFunctionSpace->meshPartition()->rankSubset();
      int computingRank = fiberNo % rankSubset->size();

      if (computingRank == rankSubset->ownRankNo())
      {
        LOG(DEBUG) << "compute (i,j)=(" << i << "," << j << "), computingRank " << computingRank 
          << ", fiberNo: " << fiberNo << ", fiberDataNo: " << fiberDataNo 
          << ", rankSubset->size(): " << rankSubset->size() << ", own: " << rankSubset->ownRankNo();
        nInstancesToCompute_ += fiberFunctionSpace->nDofsGlobal();

        assert(fiberDataNo < fiberData_.size());
        assert(fiberFunctionSpace);
        assert(motorUnitNo_.size() > 0);
        
        fiberData_.at(fiberDataNo).valuesLength = fiberFunctionSpace->nDofsGlobal();
        fiberData_.at(fiberDataNo).fiberNoGlobal = fiberNoGlobal;
        fiberData_.at(fiberDataNo).motorUnitNo = motorUnitNo_[fiberNoGlobal % motorUnitNo_.size()];
        
        // copy settings
        fiberData_.at(fiberDataNo).setSpecificStatesCallFrequency = cellmlAdapter.setSpecificStatesCallFrequency_;
        fiberData_.at(fiberDataNo).setSpecificStatesFrequencyJitter = cellmlAdapter.setSpecificStatesFrequencyJitter_;
        fiberData_.at(fiberDataNo).setSpecificStatesRepeatAfterFirstCall = cellmlAdapter.setSpecificStatesRepeatAfterFirstCall_;
        fiberData_.at(fiberDataNo).setSpecificStatesCallEnableBegin = cellmlAdapter.setSpecificStatesCallEnableBegin_;

        fiberData_.at(fiberDataNo).currentJitter = 0;
        fiberData_.at(fiberDataNo).jitterIndex = 0;
        fiberData_.at(fiberDataNo).lastStimulationCheckTime = fiberData_.at(fiberDataNo).setSpecificStatesCallEnableBegin - 1e-13 - 1./(fiberData_.at(fiberDataNo).setSpecificStatesCallFrequency+fiberData_.at(fiberDataNo).currentJitter);

        fiberData_.at(fiberDataNo).valuesOffset = 0;
        fiberData_.at(fiberDataNo).currentlyStimulating = false;
        if (fiberDataNo > 0)
        {
          fiberData_.at(fiberDataNo).valuesOffset = fiberData_.at(fiberDataNo-1).valuesOffset + fiberData_.at(fiberDataNo-1).valuesLength;
        }

        // increase index for fiberData_ struct
        fiberDataNo++;
      }
    }
  }

  // get the states and intermediates no.s to be transferred as output connector data
  CellmlAdapterType &cellmlAdapter = instances[0].timeStepping1().instancesLocal()[0].discretizableInTime();
  cellmlAdapter.getStatesIntermediatesForTransfer(statesForTransfer_, intermediatesForTransfer_);

  int nVcVectors = (nInstancesToCompute_ + Vc::double_v::Size - 1) / Vc::double_v::Size;

  fiberPointBuffers_.resize(nVcVectors);
  fiberPointBuffersIntermediatesForTransfer_.resize(nVcVectors);

  for (int i = 0; i < fiberPointBuffersIntermediatesForTransfer_.size(); i++)
  {
    fiberPointBuffersIntermediatesForTransfer_[i].resize(intermediatesForTransfer_.size());
  }

  LOG(DEBUG) << nInstancesToCompute_ << " instances to compute, " << nVcVectors
    << " Vc vectors, size of double_v: " << Vc::double_v::Size << ", "
    << statesForTransfer_.size()-1 << " additional states for transfer, "
    << intermediatesForTransfer_.size() << " intermediates for transfer";

  // initialize values
  for (int i = 0; i < fiberPointBuffers_.size(); i++)
  {
    initializeStates(fiberPointBuffers_[i].states);
  }

  // initialize field variable names
  for (int i = 0; i < instances.size(); i++)
  {
    std::vector<TimeSteppingScheme::Heun<CellmlAdapterType>> &innerInstances
      = instances[i].timeStepping1().instancesLocal();  // TimeSteppingScheme::Heun<CellmlAdapter...

    for (int j = 0; j < innerInstances.size(); j++, fiberNo++)
    {
      std::shared_ptr<FiberFunctionSpace> fiberFunctionSpace = innerInstances[j].data().functionSpace();
      CellmlAdapterType &cellmlAdapter = innerInstances[j].discretizableInTime();

      // loop over further states to transfer
      int furtherDataIndex = 0;
      for (int stateIndex = 1; stateIndex < statesForTransfer_.size(); stateIndex++, furtherDataIndex++)
      {
        // get field variable
        std::shared_ptr<FieldVariable::FieldVariable<FiberFunctionSpace,1>> fieldVariableStates
          = instances[i].timeStepping2().instancesLocal()[j].getOutputConnectorData()->variable1[stateIndex].values;

        std::string name = cellmlAdapter.data().states()->componentName(statesForTransfer_[stateIndex]);

        fieldVariableStates->setName(name);
      }

      // loop over intermediates to transfer
      for (int intermediateIndex = 0; intermediateIndex < intermediatesForTransfer_.size(); intermediateIndex++, furtherDataIndex++)
      {
        // get field variable
        std::shared_ptr<FieldVariable::FieldVariable<FiberFunctionSpace,1>> fieldVariableIntermediates
          = instances[i].timeStepping2().instancesLocal()[j].getOutputConnectorData()->variable2[intermediateIndex].values;

        std::string name = cellmlAdapter.data().intermediates()->componentName(intermediatesForTransfer_[intermediateIndex]);

        fieldVariableIntermediates->setName(name);
      }
    }
  }

}

template<int nStates, int nIntermediates>
void FastMonodomainSolverBase<nStates,nIntermediates>::
run()
{
  initialize();
  advanceTimeSpan();
}

template<int nStates, int nIntermediates>
void FastMonodomainSolverBase<nStates,nIntermediates>::
advanceTimeSpan()
{
  LOG(TRACE) << "FastMonodomainSolver::advanceTimeSpan";

  // loop over fibers and communicate element lengths and initial values to the ranks that participate in computing
  fetchFiberData();

  //Control::PerformanceMeasurement::startFlops();

  // do computation of own fibers, HH is hardcoded, stimulation from parsed MU and firing_times files
  computeMonodomain();

  //Control::PerformanceMeasurement::endFlops();

  // loop over fibers and communicate resulting values back
  updateFiberData();

  //nestedSolvers_.advanceTimeSpan();

  // call output writer of diffusion
  std::vector<typename NestedSolversType::TimeSteppingSchemeType> &instances = nestedSolvers_.instancesLocal();

  for (int i = 0; i < instances.size(); i++)
  {
    // call write output of MultipleInstances, callCountIncrement is the number of times the output writer would have been called without FastMonodomainSolver
    instances[i].timeStepping2().writeOutput(0, currentTime_, nTimeStepsSplitting_);
  }

  // call own output writers, write current 0D output values, not yet implemented
  //this->outputWriterManager_.writeOutput(*this->data_, 0, currentTime_);
}

//! get element lengths and vmValues from the other ranks
template<int nStates, int nIntermediates>
void FastMonodomainSolverBase<nStates,nIntermediates>::
fetchFiberData()
{
  VLOG(1) << "fetchFiberData";
  std::vector<typename NestedSolversType::TimeSteppingSchemeType> &instances = nestedSolvers_.instancesLocal();

  // loop over fibers and communicate element lengths and initial values to the ranks that participate in computing

  int fiberNo = 0;
  int fiberDataNo = 0;
  for (int i = 0; i < instances.size(); i++)
  {
    std::vector<TimeSteppingScheme::Heun<CellmlAdapterType>> &innerInstances
      = instances[i].timeStepping1().instancesLocal();  // TimeSteppingScheme::Heun<CellmlAdapter...

    for (int j = 0; j < innerInstances.size(); j++, fiberNo++)
    {
      std::shared_ptr<FiberFunctionSpace> fiberFunctionSpace = innerInstances[j].data().functionSpace();
      LOG(DEBUG) << "(" << i << "," << j << ") functionSpace " << fiberFunctionSpace->meshName();

      // communicate element lengths
      std::vector<double> localLengths(fiberFunctionSpace->nElementsLocal());

      // loop over local elements and compute element lengths
      for (element_no_t elementNoLocal = 0; elementNoLocal < fiberFunctionSpace->nElementsLocal(); elementNoLocal++)
      {
        std::array<Vec3, FiberFunctionSpace::nDofsPerElement()> geometryElementValues;
        fiberFunctionSpace->geometryField().getElementValues(elementNoLocal, geometryElementValues);
        double elementLength = MathUtility::distance<3>(geometryElementValues[0], geometryElementValues[1]);
        localLengths[elementNoLocal] = elementLength;
      }

      std::shared_ptr<Partition::RankSubset> rankSubset = fiberFunctionSpace->meshPartition()->rankSubset();
      MPI_Comm mpiCommunicator = rankSubset->mpiCommunicator();
      int computingRank = fiberNo % rankSubset->size();

      std::vector<int> nElementsOnRanks(rankSubset->size());
      std::vector<int> nDofsOnRanks(rankSubset->size());
      std::vector<int> offsetsOnRanks(rankSubset->size());

      double *elementLengthsReceiveBuffer = nullptr;
      double *vmValuesReceiveBuffer = nullptr;

      if (computingRank == rankSubset->ownRankNo())
      {
        fiberData_[fiberDataNo].elementLengths.resize(fiberFunctionSpace->nElementsGlobal());
        fiberData_[fiberDataNo].vmValues.resize(fiberFunctionSpace->nDofsGlobal());

        // resize buffer of further data that will be transferred back in updateFiberData()
        int nStatesAndIntermediatesValues = statesForTransfer_.size() + intermediatesForTransfer_.size() - 1;
        fiberData_[fiberDataNo].furtherStatesAndIntermediatesValues.resize(fiberFunctionSpace->nDofsGlobal() * nStatesAndIntermediatesValues);
        
        elementLengthsReceiveBuffer = fiberData_[fiberDataNo].elementLengths.data();
        vmValuesReceiveBuffer = fiberData_[fiberDataNo].vmValues.data();
      }

      for (int rankNo = 0; rankNo < rankSubset->size(); rankNo++)
      {
        nElementsOnRanks[rankNo] = fiberFunctionSpace->meshPartition()->nNodesLocalWithGhosts(0, rankNo) - 1;
        offsetsOnRanks[rankNo] = fiberFunctionSpace->meshPartition()->beginNodeGlobalNatural(0, rankNo);
        nDofsOnRanks[rankNo] = fiberFunctionSpace->meshPartition()->nNodesLocalWithoutGhosts(0, rankNo);
      }

      // int MPI_Gatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
      //          void *recvbuf, const int *recvcounts, const int *displs,
      //          MPI_Datatype recvtype, int root, MPI_Comm comm)
      //
      VLOG(1) << "Gatherv of element lengths to rank " << computingRank << ", values " << localLengths << ", sizes: " << nElementsOnRanks << ", offsets: " << offsetsOnRanks;

      MPI_Gatherv(localLengths.data(), fiberFunctionSpace->nElementsLocal(), MPI_DOUBLE,
                  elementLengthsReceiveBuffer, nElementsOnRanks.data(), offsetsOnRanks.data(),
                  MPI_DOUBLE, computingRank, mpiCommunicator);

      // communicate Vm values
      std::vector<double> vmValuesLocal;
      innerInstances[j].data().solution()->getValuesWithoutGhosts(0, vmValuesLocal);

      VLOG(1) << "Gatherv of values to rank " << computingRank << ", sizes: " << nDofsOnRanks << ", offsets: " << offsetsOnRanks << ", local values " << vmValuesLocal;

      MPI_Gatherv(vmValuesLocal.data(), fiberFunctionSpace->nDofsLocalWithoutGhosts(), MPI_DOUBLE,
                  vmValuesReceiveBuffer, nDofsOnRanks.data(), offsetsOnRanks.data(),
                  MPI_DOUBLE, computingRank, mpiCommunicator);

      // increase index for fiberData_ struct
      if (computingRank == rankSubset->ownRankNo())
        fiberDataNo++;
    }
  }

  // copy Vm values to compute buffers
  for (int fiberDataNo = 0; fiberDataNo < fiberData_.size(); fiberDataNo++)
  {
    int nValues = fiberData_[fiberDataNo].vmValues.size();

    for (int valueNo = 0; valueNo < nValues; valueNo++)
    {
      global_no_t valueIndexAllFibers = fiberData_[fiberDataNo].valuesOffset + valueNo;

      global_no_t pointBuffersNo = valueIndexAllFibers / Vc::double_v::Size;
      int entryNo = valueIndexAllFibers % Vc::double_v::Size;

      fiberPointBuffers_[pointBuffersNo].states[0][entryNo] = fiberData_[fiberDataNo].vmValues[valueNo];
    }
  }
}

//! send vmValues data from fiberData_ back to the fibers where it belongs to and set in the respective field variable
template<int nStates, int nIntermediates>
void FastMonodomainSolverBase<nStates,nIntermediates>::
updateFiberData()
{
  // copy Vm and other states/intermediates from compute buffers to fiberData_
  for (int fiberDataNo = 0; fiberDataNo < fiberData_.size(); fiberDataNo++)
  {
    int nValues = fiberData_[fiberDataNo].vmValues.size();

    for (int valueNo = 0; valueNo < nValues; valueNo++)
    {
      global_no_t valueIndexAllFibers = fiberData_[fiberDataNo].valuesOffset + valueNo;

      global_no_t pointBuffersNo = valueIndexAllFibers / Vc::double_v::Size;
      int entryNo = valueIndexAllFibers % Vc::double_v::Size;

      assert(statesForTransfer_.size() > 0);
      const int stateToTransfer = statesForTransfer_[0];  // transfer the first state value
      fiberData_[fiberDataNo].vmValues[valueNo] = fiberPointBuffers_[pointBuffersNo].states[stateToTransfer][entryNo];

      // loop over further states to transfer
      int furtherDataIndex = 0;
      for (int i = 1; i < statesForTransfer_.size(); i++, furtherDataIndex++)
      {
        const int stateToTransfer = statesForTransfer_[i];

        fiberData_[fiberDataNo].furtherStatesAndIntermediatesValues[furtherDataIndex*nValues + valueNo]
          = fiberPointBuffers_[pointBuffersNo].states[stateToTransfer][entryNo];
      }

      // loop over intermediates to transfer
      for (int i = 0; i < intermediatesForTransfer_.size(); i++, furtherDataIndex++)
      {
        fiberData_[fiberDataNo].furtherStatesAndIntermediatesValues[furtherDataIndex*nValues + valueNo]
          = fiberPointBuffersIntermediatesForTransfer_[pointBuffersNo][i][entryNo];
      }
    }
    LOG(DEBUG) << "states and intermediates for transfer at fiberDataNo=" << fiberDataNo << ": " << fiberData_[fiberDataNo].furtherStatesAndIntermediatesValues;
    LOG(DEBUG) << "size: " << fiberData_[fiberDataNo].furtherStatesAndIntermediatesValues.size() << ", nValues: " << nValues;
  }

  LOG(TRACE) << "updateFiberData";
  std::vector<typename NestedSolversType::TimeSteppingSchemeType> &instances = nestedSolvers_.instancesLocal();

  // loop over fibers and communicate element lengths and initial values to the ranks that participate in computing

  int fiberNo = 0;
  int fiberDataNo = 0;
  for (int i = 0; i < instances.size(); i++)
  {
    std::vector<TimeSteppingScheme::Heun<CellmlAdapterType>> &innerInstances
      = instances[i].timeStepping1().instancesLocal();  // TimeSteppingScheme::Heun<CellmlAdapter...

    for (int j = 0; j < innerInstances.size(); j++, fiberNo++)
    {
      std::shared_ptr<FiberFunctionSpace> fiberFunctionSpace = innerInstances[j].data().functionSpace();
      //CellmlAdapterType &cellmlAdapter = innerInstances[j].discretizableInTime();

      // prepare helper variables for Scatterv
      std::shared_ptr<Partition::RankSubset> rankSubset = fiberFunctionSpace->meshPartition()->rankSubset();
      MPI_Comm mpiCommunicator = rankSubset->mpiCommunicator();
      int computingRank = fiberNo % rankSubset->size();

      std::vector<int> nDofsOnRanks(rankSubset->size());
      std::vector<int> offsetsOnRanks(rankSubset->size());

      for (int rankNo = 0; rankNo < rankSubset->size(); rankNo++)
      {
        offsetsOnRanks[rankNo] = fiberFunctionSpace->meshPartition()->beginNodeGlobalNatural(0, rankNo);
        nDofsOnRanks[rankNo] = fiberFunctionSpace->meshPartition()->nNodesLocalWithoutGhosts(0, rankNo);
      }

      double *sendBufferVmValues = nullptr;
      if (computingRank == rankSubset->ownRankNo())
      {
        sendBufferVmValues = fiberData_[fiberDataNo].vmValues.data();
      }

      //  int MPI_Scatterv(const void *sendbuf, const int *sendcounts, const int *displs, MPI_Datatype sendtype,
      //                  void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
      // communicate Vm values
      std::vector<double> vmValuesLocal(fiberFunctionSpace->nDofsLocalWithoutGhosts());
      MPI_Scatterv(sendBufferVmValues, nDofsOnRanks.data(), offsetsOnRanks.data(), MPI_DOUBLE,
                   vmValuesLocal.data(), fiberFunctionSpace->nDofsLocalWithoutGhosts(), MPI_DOUBLE,
                   computingRank, mpiCommunicator);

      LOG(DEBUG) << "Scatterv from rank " << computingRank << ", sizes: " << nDofsOnRanks << ", offsets: " << offsetsOnRanks << ", received local values " << vmValuesLocal;

      // store Vm values in CellmlAdapter and diffusion FiniteElementMethod
      LOG(DEBUG) << "fiber " << fiberDataNo << ", set values " << vmValuesLocal;
      innerInstances[j].data().solution()->setValuesWithoutGhosts(0, vmValuesLocal);
      instances[i].timeStepping2().instancesLocal()[j].data().solution()->setValuesWithoutGhosts(0, vmValuesLocal);

      // ----------------------
      // communicate further states and intermediates that are selected by the options "statesForTransfer" and "intermediatesForTransfer"

      std::vector<int> nValuesOnRanks(rankSubset->size());
      int nStatesAndIntermediatesValues = statesForTransfer_.size() + intermediatesForTransfer_.size() - 1;

      for (int rankNo = 0; rankNo < rankSubset->size(); rankNo++)
      {
        offsetsOnRanks[rankNo] = fiberFunctionSpace->meshPartition()->beginNodeGlobalNatural(0, rankNo) * nStatesAndIntermediatesValues;
        nValuesOnRanks[rankNo] = fiberFunctionSpace->meshPartition()->nNodesLocalWithoutGhosts(0, rankNo) * nStatesAndIntermediatesValues;
      }

      double *sendBuffer = nullptr;
      if (computingRank == rankSubset->ownRankNo())
      {
        sendBuffer = fiberData_[fiberDataNo].furtherStatesAndIntermediatesValues.data();
      }

      std::vector<double> valuesLocal(fiberFunctionSpace->nDofsLocalWithoutGhosts() * nStatesAndIntermediatesValues);
      MPI_Scatterv(sendBuffer, nValuesOnRanks.data(), offsetsOnRanks.data(), MPI_DOUBLE,
                   valuesLocal.data(), fiberFunctionSpace->nDofsLocalWithoutGhosts() * nStatesAndIntermediatesValues, MPI_DOUBLE,
                   computingRank, mpiCommunicator);

      LOG(DEBUG) << "Scatterv furtherStatesAndIntermediatesValues from rank " << computingRank << ", sizes: " << nValuesOnRanks << ", offsets: " << offsetsOnRanks
        << ", sendBuffer: " << sendBuffer << ", received local values: " << valuesLocal;

      // store received states and intermediates values in diffusion outputConnectorData
      // loop over further states to transfer
      int furtherDataIndex = 0;
      for (int stateIndex = 1; stateIndex < statesForTransfer_.size(); stateIndex++, furtherDataIndex++)
      {
        // store in diffusion
        std::shared_ptr<FieldVariable::FieldVariable<FiberFunctionSpace,1>> fieldVariableStates
          = instances[i].timeStepping2().instancesLocal()[j].getOutputConnectorData()->variable1[stateIndex].values;

        int nValues = fiberFunctionSpace->nDofsLocalWithoutGhosts();
        double *values = valuesLocal.data() + furtherDataIndex * nValues;

        // int componentNo, int nValues, const dof_no_t *dofNosLocal, const double *values
        fieldVariableStates->setValues(0, nValues, fiberFunctionSpace->meshPartition()->dofNosLocal().data(), values);

        // store in cellmlAdapter
        std::shared_ptr<FieldVariable::FieldVariable<FiberFunctionSpace,nStates>> fieldVariableStatesCellML
          = instances[i].timeStepping1().instancesLocal()[j].getOutputConnectorData()->variable1[stateIndex].values;

        const int componentNo = statesForTransfer_[stateIndex];

        // int componentNo, int nValues, const dof_no_t *dofNosLocal, const double *values
        fieldVariableStatesCellML->setValues(componentNo, nValues, fiberFunctionSpace->meshPartition()->dofNosLocal().data(), values);

        LOG(DEBUG) << "store " << nValues << " values for additional state " << statesForTransfer_[stateIndex];
      }

      // loop over intermediates to transfer
      for (int intermediateIndex = 0; intermediateIndex < intermediatesForTransfer_.size(); intermediateIndex++, furtherDataIndex++)
      {
        // store in diffusion
        std::shared_ptr<FieldVariable::FieldVariable<FiberFunctionSpace,1>> fieldVariableIntermediates
          = instances[i].timeStepping2().instancesLocal()[j].getOutputConnectorData()->variable2[intermediateIndex].values;

        int nValues = fiberFunctionSpace->nDofsLocalWithoutGhosts();
        double *values = valuesLocal.data() + furtherDataIndex * nValues;

        // int componentNo, int nValues, const dof_no_t *dofNosLocal, const double *values
        fieldVariableIntermediates->setValues(0, nValues, fiberFunctionSpace->meshPartition()->dofNosLocal().data(), values);

        // store in CellmlAdapter
        std::shared_ptr<FieldVariable::FieldVariable<FiberFunctionSpace,1>> fieldVariableIntermediatesCellML
          = instances[i].timeStepping1().instancesLocal()[j].getOutputConnectorData()->variable2[intermediateIndex].values;

        //const int componentNo = intermediatesForTransfer_[intermediateIndex];

        // int componentNo, int nValues, const dof_no_t *dofNosLocal, const double *values
        fieldVariableIntermediatesCellML->setValues(0, nValues, fiberFunctionSpace->meshPartition()->dofNosLocal().data(), values);

        LOG(DEBUG) << "store " << nValues << " values for intermediate " << intermediatesForTransfer_[intermediateIndex];
        LOG(DEBUG) << *fieldVariableIntermediates;
      }

      // increase index for fiberData_ struct
      if (computingRank == rankSubset->ownRankNo())
        fiberDataNo++;
    }
  }
}

template<int nStates, int nIntermediates>
void FastMonodomainSolverBase<nStates,nIntermediates>::
computeMonodomain()
{
  LOG(TRACE) << "computeMonodomain";

  // initialize data vector
  // array of vectorized struct

  // fetch timestep widths and total time span
  std::vector<typename NestedSolversType::TimeSteppingSchemeType> &instances = nestedSolvers_.instancesLocal();

  TimeSteppingScheme::Heun<CellmlAdapterType> &heun = instances[0].timeStepping1().instancesLocal()[0];
  durationLogKey0D_ = heun.durationLogKey();

  ImplicitEuler &implicitEuler = instances[0].timeStepping2().instancesLocal()[0];
  durationLogKey1D_ = implicitEuler.durationLogKey();
  double prefactor = implicitEuler.discretizableInTime().data().context().getPythonConfig().getOptionDouble("prefactor", 1.0);

  LOG(DEBUG) << "durationLogKeys: " << durationLogKey0D_ << "," << durationLogKey1D_;

  double startTime = instances[0].startTime();
  double timeStepWidthSplitting = instances[0].timeStepWidth();
  nTimeStepsSplitting_ = instances[0].numberTimeSteps();

  heun.setTimeSpan(startTime, startTime + 0.5 * timeStepWidthSplitting);
  double dt0D = heun.timeStepWidth();
  int nTimeSteps0D = heun.numberTimeSteps();

  implicitEuler.setTimeSpan(startTime, startTime + timeStepWidthSplitting);
  double dt1D = implicitEuler.timeStepWidth();
  int nTimeSteps1D = implicitEuler.numberTimeSteps();

  LOG(DEBUG) << "prefactor: " << prefactor << ", dtSplitting: " << timeStepWidthSplitting << ", n steps: " << nTimeStepsSplitting_;
  LOG(DEBUG) << "dt0D: " << dt0D << ", n steps: " << nTimeSteps0D << ", dt1D: " << dt1D << ", n steps: " << nTimeSteps1D;

  // picture for strang splitting:
  // ===== t ==>
  // -1->
  //  /
  // ---2--->
  //      /
  //     -1->
  //        |
  //        2

  // loop over splitting time steps
  for (int timeStepNo = 0; timeStepNo < nTimeStepsSplitting_; timeStepNo++)
  {
    // perform Strang splitting
    double currentTime = startTime + timeStepNo * timeStepWidthSplitting;

    LOG(DEBUG) << "splitting " << timeStepNo << "/" << nTimeStepsSplitting_ << ", t: " << currentTime;

    // compute midTime once per step to reuse it. [currentTime, midTime=currentTime+0.5*timeStepWidth, currentTime+timeStepWidth]
    double midTime = currentTime + 0.5 * timeStepWidthSplitting;
    bool storeIntermediatesForTransfer = timeStepNo == nTimeStepsSplitting_-1;   // after the last timestep, store the intermediates for transfer

    // perform splitting
    compute0D(currentTime, dt0D, nTimeSteps0D, false);
    compute1D(currentTime, dt1D, nTimeSteps1D, prefactor);
    compute0D(midTime,     dt0D, nTimeSteps0D, storeIntermediatesForTransfer);
  }

  currentTime_ = instances[0].endTime();
}

template<int nStates, int nIntermediates>
void FastMonodomainSolverBase<nStates,nIntermediates>::
compute0D(double startTime, double timeStepWidth, int nTimeSteps, bool storeIntermediatesForTransfer)
{
  Control::PerformanceMeasurement::start(durationLogKey0D_);
  LOG(DEBUG) << "compute0D(" << startTime << "), " << nTimeSteps << " time steps";

  using Vc::double_v;

  // Heun scheme:
  // y* = y_n + dt*rhs(y_n)
  // y_n+1 = y_n + 0.5*[rhs(y_n) + rhs(y*)]

  // loop over buffers
  for (global_no_t pointBuffersNo = 0; pointBuffersNo < fiberPointBuffers_.size(); pointBuffersNo++)
  {
    int fiberDataNo = pointBuffersNo * Vc::double_v::Size / fiberData_[0].valuesLength;
    int indexInFiber = pointBuffersNo * Vc::double_v::Size - fiberData_[fiberDataNo].valuesOffset;
    int fiberCenterIndex = fiberData_[fiberDataNo].valuesLength / 2;
    bool currentPointIsInCenter = ((fiberCenterIndex - indexInFiber) < Vc::double_v::Size);
    VLOG(3) << "currentPointIsInCenter: " << currentPointIsInCenter << ", pointBuffersNo: " << pointBuffersNo << ", fiberDataNo: " << fiberDataNo << ", indexInFiber:" << indexInFiber << ", fiberCenterIndex: " << fiberCenterIndex << ", " << (indexInFiber - fiberCenterIndex) << " < " << Vc::double_v::Size;

    const int &motorUnitNo = fiberData_[fiberDataNo].motorUnitNo;
    VLOG(3) << "pointBuffersNo: " << pointBuffersNo << ", fiberDataNo: " << fiberDataNo << ", indexInFiber: " << indexInFiber << ", motorUnitNo: " << motorUnitNo;

    /*
    double_v state0 = fiberPointBuffers_[pointBuffersNo].states[0];
    double_v state1 = fiberPointBuffers_[pointBuffersNo].states[1];
    double_v state2 = fiberPointBuffers_[pointBuffersNo].states[2];
    double_v state3 = fiberPointBuffers_[pointBuffersNo].states[3];

    VLOG(3) << "  states [" << state0 << "," << state1 << "," << state2 << "," << state3 << "]";
*/
    // loop over timesteps
    for (int timeStepNo = 0; timeStepNo < nTimeSteps; timeStepNo++)
    {
      // determine if fiber gets stimulated
      double currentTime = startTime + timeStepNo * timeStepWidth;

      // get time from with testing for stimulation is enabled
      const double lastStimulationCheckTime = fiberData_[fiberDataNo].lastStimulationCheckTime;

      const double setSpecificStatesCallFrequency = fiberData_[fiberDataNo].setSpecificStatesCallFrequency;
      const double setSpecificStatesRepeatAfterFirstCall = fiberData_[fiberDataNo].setSpecificStatesRepeatAfterFirstCall;
      const double setSpecificStatesCallEnableBegin = fiberData_[fiberDataNo].setSpecificStatesCallEnableBegin;

      std::vector<double> &setSpecificStatesFrequencyJitter = fiberData_[fiberDataNo].setSpecificStatesFrequencyJitter;
      int &jitterIndex = fiberData_[fiberDataNo].jitterIndex;
      double &currentJitter = fiberData_[fiberDataNo].currentJitter;

      // check if time has come to call setSpecificStates
      bool checkStimulation = false;

      VLOG(1) << "currentTime: " << currentTime << ", lastStimulationCheckTime: " << lastStimulationCheckTime << ", next time point: " << lastStimulationCheckTime + 1./(setSpecificStatesCallFrequency+currentJitter);
      VLOG(1) << "setSpecificStatesCallFrequency: " << setSpecificStatesCallFrequency << ", currentJitter: " << currentJitter << ", setSpecificStatesCallEnableBegin: " << setSpecificStatesCallEnableBegin;

      if (currentTime >= lastStimulationCheckTime + 1./(setSpecificStatesCallFrequency+currentJitter)
          && currentTime >= setSpecificStatesCallEnableBegin-1e-13)
      {
        VLOG(1) << "-> checkStimulation";
        checkStimulation = true;

        VLOG(1) << "check if stimulation is over: duration already: " << currentTime - (lastStimulationCheckTime + 1./(setSpecificStatesCallFrequency+currentJitter))
          << ", setSpecificStatesRepeatAfterFirstCall: " << setSpecificStatesRepeatAfterFirstCall;

        // if current stimulation is over
        if (currentTime - (lastStimulationCheckTime + 1./(setSpecificStatesCallFrequency+currentJitter)) > setSpecificStatesRepeatAfterFirstCall)
        {
          // advance time of last call to specificStates
          LOG(DEBUG) << " old lastStimulationCheckTime: " << fiberData_[fiberDataNo].lastStimulationCheckTime << ", currentJitter: " << currentJitter << ", add " << 1./(setSpecificStatesCallFrequency+currentJitter);
          fiberData_[fiberDataNo].lastStimulationCheckTime += 1./(setSpecificStatesCallFrequency+currentJitter);

          LOG(DEBUG) << " new lastStimulationCheckTime: " << fiberData_[fiberDataNo].lastStimulationCheckTime;

          // compute new jitter value
          double jitterFactor = 0.0;
          if (setSpecificStatesFrequencyJitter.size() > 0)
            jitterFactor = setSpecificStatesFrequencyJitter[jitterIndex % setSpecificStatesFrequencyJitter.size()];
          currentJitter = jitterFactor * setSpecificStatesCallFrequency;
          LOG(DEBUG) << " jitterIndex: " << jitterIndex << ", new jitterFactor: " << jitterFactor << ", currentJitter: " << currentJitter;
          jitterIndex++;

          checkStimulation = false;
        }
      }

      // instead of calling setSpecificStates, directly determine whether to stimulate from the firingEvents file
      int firingEventsIndex = round(currentTime * setSpecificStatesCallFrequency);

      bool stimulate =
        checkStimulation
        && firingEvents_[firingEventsIndex % firingEvents_.size()][motorUnitNo % firingEvents_[firingEventsIndex % firingEvents_.size()].size()];

      if (checkStimulation)
      {
        VLOG(1) << "setSpecificStatesCallFrequency: " << setSpecificStatesCallFrequency << ", firingEventsIndex: " << firingEventsIndex << ", fires: "
          << firingEvents_[firingEventsIndex % firingEvents_.size()][motorUnitNo % firingEvents_[firingEventsIndex % firingEvents_.size()].size()];
        VLOG(1) << "currentPointIsInCenter: " << currentPointIsInCenter;
      }

      if (stimulate && currentPointIsInCenter)
      {
        // if this is the first point in time of the current stimulation, log stimulation time
        if (!fiberData_[fiberDataNo].currentlyStimulating)
        {
          fiberData_[fiberDataNo].currentlyStimulating = true;
          Control::StimulationLogging::logStimulationBegin(currentTime, fiberData_[fiberDataNo].motorUnitNo, fiberData_[fiberDataNo].fiberNoGlobal);
        }

        LOG(DEBUG) << "stimulate fiber " << fiberData_[fiberDataNo].fiberNoGlobal << ", MU " << motorUnitNo << " at t=" << currentTime;
        LOG(DEBUG) << "  pointBuffersNo: " << pointBuffersNo << ", indexInFiber: " << indexInFiber << ", fiberCenterIndex: " << fiberCenterIndex;
        LOG(DEBUG) << "  motorUnitNo: " << motorUnitNo << " (" << motorUnitNo % firingEvents_[firingEventsIndex % firingEvents_.size()].size() << ")";
        LOG(DEBUG) << "  firing events index: " << firingEventsIndex << " (" << firingEventsIndex % firingEvents_.size() << ")";
        LOG(DEBUG) << "  setSpecificStatesCallEnableBegin: " << setSpecificStatesCallEnableBegin << ", lastStimulationCheckTime: " << lastStimulationCheckTime
          << ", stimulation already for " << + 1./(setSpecificStatesCallFrequency+currentJitter);
      }

      if (!stimulate)
      {
        fiberData_[fiberDataNo].currentlyStimulating = false;
      }

      if (stimulate && currentPointIsInCenter)
        LOG(INFO) << "t: " << currentTime << ", stimulate fiber " << fiberData_[fiberDataNo].fiberNoGlobal << ", MU " << motorUnitNo;

      const bool argumentStimulate = stimulate && currentPointIsInCenter;
      const bool argumentStoreIntermediates = storeIntermediatesForTransfer && timeStepNo == nTimeSteps-1;

      /*LOG(DEBUG) << "storeIntermediatesForTransfer: " << storeIntermediatesForTransfer << ", timeStepNo " << timeStepNo << ", nTimeSteps: " << nTimeSteps
        << ", argumentStoreIntermediates: " << argumentStoreIntermediates;*/

      compute0DInstance(fiberPointBuffers_[pointBuffersNo].states, currentTime, timeStepWidth, argumentStimulate,
                        argumentStoreIntermediates, fiberPointBuffersIntermediatesForTransfer_[pointBuffersNo]);

      // perform one step of the heun scheme

    }


    //VLOG(3) << "-> index " << pointBuffersNo << ", states [" << state0 << "," << state1 << "," << state2 << "," << state3 << "]";
  }

  VLOG(1) << "nFiberPointBuffers: " << fiberPointBuffers_.size();
  Control::PerformanceMeasurement::stop(durationLogKey0D_);
}

void compute0DInstance();

template<int nStates, int nIntermediates>
void FastMonodomainSolverBase<nStates,nIntermediates>::
compute1D(double startTime, double timeStepWidth, int nTimeSteps, double prefactor)
{
  Control::PerformanceMeasurement::start(durationLogKey1D_);

  LOG(DEBUG) << "compute1D(" << startTime << ")";

  // perform implicit euler step
  // (K - 1/dt*M) u^{n+1} = -1/dt*M u^{n})

  // stencil K: 1/h*[_-1_  1  ]*prefactor
  // stencil M:   h*[_1/3_ 1/6]

  const double dt = timeStepWidth;

  // loop over fibers
  for (int fiberDataNo = 0; fiberDataNo < fiberData_.size(); fiberDataNo++)
  {
    int nValues = fiberData_[fiberDataNo].vmValues.size();

#ifndef NDEBUG
    VLOG(1) << "fiber " << fiberDataNo << "/" << fiberData_.size() << ", valuesOffset: " << fiberData_[fiberDataNo].valuesOffset
      << ", has " << nValues << " values: ";

    std::stringstream s, s2;
    for (int valueNo = 0; valueNo < nValues; valueNo++)
    {
      global_no_t valuesIndexAllFibers = fiberData_[fiberDataNo].valuesOffset + valueNo;
      global_no_t pointBuffersNo = valuesIndexAllFibers / Vc::double_v::Size;
      int entryNo = valuesIndexAllFibers % Vc::double_v::Size;
      double u = fiberPointBuffers_[pointBuffersNo].states[0][entryNo];
      if (valueNo != 0)
      {
        s << ", ";
        s2 << ",";
      }
      s << u;
      s2 << "[" << pointBuffersNo << "][" << entryNo << "](" << valuesIndexAllFibers << ")  ";
    }
    VLOG(1) << s.str();
    //LOG(DEBUG) << "indices: " << s2.str();
#endif
    // [ b c     ] [x]   [d]
    // [ a b c   ] [x] = [d]
    // [   a b c ] [x]   [d]
    // [     a b ] [x]   [d]

    // Thomas algorithm
    // forward substitution
    // c'_0 = c_0 / b_0
    // c'_i = c_i / (b_i - c'_{i-1}*a_i)

    // d'_0 = d_0 / b_0
    // d'_i = (d_i - d'_{i-1}*a_i) / (b_i - c'_{i-1}*a_i)

    // backward substitution
    // x_n = d'_n
    // x_i = d'_i - c'_i * x_{i+1}

    // helper buffers c', d'
    static std::vector<double> cIntermediate(nValues-1);
    static std::vector<double> dIntermediate(nValues);

    // perform forward substitution
    // loop over entries / rows of matrices
    for (int valueNo = 0; valueNo < nValues; valueNo++)
    {
      double a = 0;
      double b = 0;
      double c = 0;
      double d = 0;

      double u_previous = 0;
      double u_center = 0;
      double u_next = 0;

      global_no_t valuesIndexAllFibers = fiberData_[fiberDataNo].valuesOffset + valueNo;
      global_no_t pointBuffersNo = valuesIndexAllFibers / Vc::double_v::Size;
      int entryNo = valuesIndexAllFibers % Vc::double_v::Size;
      u_center = fiberPointBuffers_[pointBuffersNo].states[0][entryNo];

      // contribution from left element
      if (valueNo > 0)
      {
        global_no_t valuesIndexAllFibers = fiberData_[fiberDataNo].valuesOffset + valueNo - 1;
        global_no_t pointBuffersNo = valuesIndexAllFibers / Vc::double_v::Size;
        int entryNo = valuesIndexAllFibers % Vc::double_v::Size;
        u_previous = fiberPointBuffers_[pointBuffersNo].states[0][entryNo];

        // stencil K: 1/h*[1   _-1_ ]*prefactor
        // stencil M:   h*[1/6 _1/3_]

        double h_left = fiberData_[fiberDataNo].elementLengths[valueNo-1];
        double k_left = 1./h_left*(1) * prefactor;
        double m_left = h_left*1./6;
        a = (k_left - 1/dt*m_left);

        double k_right = 1./h_left*(-1) * prefactor;
        double m_right = h_left*1./3;
        b += (k_right - 1/dt*m_right);

        d += (-1/dt*m_left) * u_previous + (-1/dt*m_right) * u_center;
      }

      // contribution from right element
      if (valueNo < nValues-1)
      {
        global_no_t valuesIndexAllFibers = fiberData_[fiberDataNo].valuesOffset + valueNo + 1;
        global_no_t pointBuffersNo = valuesIndexAllFibers / Vc::double_v::Size;
        int entryNo = valuesIndexAllFibers % Vc::double_v::Size;
        u_next = fiberPointBuffers_[pointBuffersNo].states[0][entryNo];

        // stencil K: 1/h*[_-1_  1  ]*prefactor
        // stencil M:   h*[_1/3_ 1/6]

        double h_right = fiberData_[fiberDataNo].elementLengths[valueNo];
        double k_right = 1./h_right*(1) * prefactor;
        double m_right = h_right*1./6;
        c = (k_right - 1/dt*m_right);

        double k_left = 1./h_right*(-1) * prefactor;
        double m_left = h_right*1./3;
        b += (k_left - 1/dt*m_left);

        d += (-1/dt*m_left) * u_center + (-1/dt*m_right) * u_next;
      }

      if (valueNo == 0)
      {
        // c'_0 = c_0 / b_0
        cIntermediate[valueNo] = c / b;

        // d'_0 = d_0 / b_0
        dIntermediate[valueNo] = d / b;
      }
      else
      {
        if (valueNo != nValues-1)
        {
          // c'_i = c_i / (b_i - c'_{i-1}*a_i)
          cIntermediate[valueNo] = c / (b - cIntermediate[valueNo-1]*a);
        }

        // d'_i = (d_i - d'_{i-1}*a_i) / (b_i - c'_{i-1}*a_i)
        dIntermediate[valueNo] = (d - dIntermediate[valueNo-1]*a) / (b - cIntermediate[valueNo-1]*a);
      }

#ifndef NDEBUG
      if (valueNo < 5 || valueNo >= nValues-6)
      {
        VLOG(2) << "valueNo: " << valueNo << ", a: " << a << ", b: " << b << ", c: " << c << ", d: " << d
          << ", u: " << u_center << ", c': " << (valueNo < nValues-1? cIntermediate[valueNo] : 0) << ", d': " << dIntermediate[valueNo]
          << " = (" << d << "-" << dIntermediate[valueNo-1]*a << ")/(" << b << "-" << cIntermediate[valueNo-1]*a << ") = " << (d - dIntermediate[valueNo-1]*a) << "/" << (b - cIntermediate[valueNo-1]*a);
      }
#endif

    }

    //LOG(DEBUG) << "cIntermediate: " << cIntermediate;
    //LOG(DEBUG) << "dIntermediate: " << dIntermediate;

    // perform backward substitution
    // x_n = d'_n
    global_no_t valuesIndexAllFibers = fiberData_[fiberDataNo].valuesOffset + nValues-1;
    global_no_t pointBuffersNo = valuesIndexAllFibers / Vc::double_v::Size;
    int entryNo = valuesIndexAllFibers % Vc::double_v::Size;
    fiberPointBuffers_[pointBuffersNo].states[0][entryNo] = dIntermediate[nValues-1];
    //LOG(DEBUG) << "set entry (" << pointBuffersNo << "," << entryNo << ") to " << dIntermediate[nValues-1];

    double previousValue = dIntermediate[nValues-1];

    // loop over entries / rows of matrices
    for (int valueNo = nValues-2; valueNo >= 0; valueNo--)
    {
      // x_i = d'_i - c'_i * x_{i+1}
      global_no_t valuesIndexAllFibers = fiberData_[fiberDataNo].valuesOffset + valueNo;
      global_no_t pointBuffersNo = valuesIndexAllFibers / Vc::double_v::Size;
      int entryNo = valuesIndexAllFibers % Vc::double_v::Size;

      double resultValue = dIntermediate[valueNo] - cIntermediate[valueNo] * previousValue;
      fiberPointBuffers_[pointBuffersNo].states[0][entryNo] = resultValue;

      previousValue = resultValue;
    }

#ifndef NDEBUG
    s.str("");
    for (int valueNo = 0; valueNo < nValues; valueNo++)
    {
      global_no_t valuesIndexAllFibers = fiberData_[fiberDataNo].valuesOffset + valueNo;
      global_no_t pointBuffersNo = valuesIndexAllFibers / Vc::double_v::Size;
      int entryNo = valuesIndexAllFibers % Vc::double_v::Size;
      double u = fiberPointBuffers_[pointBuffersNo].states[0][entryNo];
      if (valueNo != 0)
        s << ", ";
      s << u;
    }
    VLOG(1) << " -> " << s.str();
#endif
  }
  Control::PerformanceMeasurement::stop(durationLogKey1D_);
}

template<int nStates, int nIntermediates>
void FastMonodomainSolverBase<nStates,nIntermediates>::
reset()
{
  nestedSolvers_.reset();
}

template<int nStates, int nIntermediates>
typename FastMonodomainSolverBase<nStates,nIntermediates>::Data &FastMonodomainSolverBase<nStates,nIntermediates>::
data()
{
  return nestedSolvers_.data();
}

template<int nStates, int nIntermediates>
void FastMonodomainSolverBase<nStates,nIntermediates>::
setTimeSpan(double startTime, double endTime)
{
  nestedSolvers_.setTimeSpan(startTime, endTime);
}

template<int nStates, int nIntermediates>
std::shared_ptr<typename FastMonodomainSolverBase<nStates,nIntermediates>::OutputConnectorDataType>
FastMonodomainSolverBase<nStates,nIntermediates>::
getOutputConnectorData()
{
  return nestedSolvers_.getOutputConnectorData();
}
