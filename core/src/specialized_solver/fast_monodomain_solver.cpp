#include "specialized_solver/fast_monodomain_solver.h"

#include "partition/rank_subset.h"

FastMonodomainSolver<Control::MultipleInstances<OperatorSplitting::Strang<Control::MultipleInstances<TimeSteppingScheme::Heun<CellmlAdapter<4, 9, FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1> > > > >, Control::MultipleInstances<TimeSteppingScheme::ImplicitEuler<SpatialDiscretization::FiniteElementMethod<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1>, Quadrature::Gauss<2>, Equation::Dynamic::IsotropicDiffusion> > > > > >::
FastMonodomainSolver(const DihuContext &context) :
  specificSettings_(context.getPythonConfig()), nestedSolvers_(context)
{
}

void FastMonodomainSolver<Control::MultipleInstances<OperatorSplitting::Strang<Control::MultipleInstances<TimeSteppingScheme::Heun<CellmlAdapter<4, 9, FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1> > > > >, Control::MultipleInstances<TimeSteppingScheme::ImplicitEuler<SpatialDiscretization::FiniteElementMethod<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1>, Quadrature::Gauss<2>, Equation::Dynamic::IsotropicDiffusion> > > > > >::
initialize()
{
  nestedSolvers_.initialize();

  // initialize motor unit numbers and firing times
  fiberDistributionFilename_ = specificSettings_.getOptionString("fiberDistributionFile", "");
  firingTimesFilename_ = specificSettings_.getOptionString("firingTimesFile", "");

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
  //LOG(DEBUG) << "firingEvents_:" << firingEvents_;
  LOG(DEBUG) << "motorUnitNo_: " << motorUnitNo_;

  if (motorUnitNo_.empty())
    LOG(FATAL) << "Could not parse motor units.";

  if (firingEvents_.empty())
    LOG(FATAL) << "Could not parse firing times.";

  // initialize data structures
  std::vector<NestedSolversType::TimeSteppingSchemeType> &instances = nestedSolvers_.instancesLocal();

  // determine number of fibers
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
      CellmlAdapterType &cellmlAdapter = innerInstances[j].discretizableInTime();
      setSpecificStatesCallFrequency_ = cellmlAdapter.setSpecificStatesCallFrequency_;
      setSpecificStatesRepeatAfterFirstCall_ = cellmlAdapter.setSpecificStatesRepeatAfterFirstCall_;

      std::shared_ptr<FiberFunctionSpace> fiberFunctionSpace = innerInstances[j].data().functionSpace();
      std::shared_ptr<Partition::RankSubset> rankSubset = fiberFunctionSpace->meshPartition()->rankSubset();
      int computingRank = fiberNo % rankSubset->size();

      if (computingRank == rankSubset->ownRankNo())
      {
        nFibersToCompute_++;
      }
    }
  }

  fiberData_.resize(nFibersToCompute_);
  LOG(DEBUG) << "nFibers: " << nFibers << ", nFibersToCompute_: " << nFibersToCompute_;

  // determine total number of Hodgkin-Huxley instances to compute on this rank
  nInstancesToCompute_ = 0;
  int fiberDataNo = 0;
  for (int i = 0; i < instances.size(); i++)
  {
    std::vector<TimeSteppingScheme::Heun<CellmlAdapterType>> &innerInstances
      = instances[i].timeStepping1().instancesLocal();  // TimeSteppingScheme::Heun<CellmlAdapter...

    for (int j = 0; j < innerInstances.size(); j++, fiberNo++)
    {
      CellmlAdapterType &cellmlAdapter = innerInstances[j].discretizableInTime();
      int fiberNo = PythonUtility::convertFromPython<int>::get(cellmlAdapter.pySetFunctionAdditionalParameter_);

      std::shared_ptr<FiberFunctionSpace> fiberFunctionSpace = innerInstances[j].data().functionSpace();

      std::shared_ptr<Partition::RankSubset> rankSubset = fiberFunctionSpace->meshPartition()->rankSubset();
      int computingRank = fiberNo % rankSubset->size();

      if (computingRank == rankSubset->ownRankNo())
      {
        nInstancesToCompute_ += fiberFunctionSpace->nDofsGlobal();
        fiberData_[fiberDataNo].valuesLength = fiberFunctionSpace->nDofsGlobal();
        fiberData_[fiberDataNo].fiberNo = fiberNo;
        fiberData_[fiberDataNo].motorUnitNo = motorUnitNo_[fiberNo % motorUnitNo_.size()];

        fiberData_[fiberDataNo].valuesOffset = 0;
        if (fiberDataNo > 0)
        {
          fiberData_[fiberDataNo].valuesOffset = fiberData_[fiberDataNo-1].valuesOffset + fiberData_[fiberDataNo].valuesLength;
        }

        // increase index for fiberData_ struct
        fiberDataNo++;
      }
    }
  }

  int nVcVectors = (nInstancesToCompute_ + Vc::double_v::Size - 1) / Vc::double_v::Size;

  fiberPointBuffers_.resize(nVcVectors);
  LOG(DEBUG) << nInstancesToCompute_ << " instances to compute, " << nVcVectors << " Vc vectors, size of double_v: " << Vc::double_v::Size;

  // initialize values
  for (int i = 0; i < fiberPointBuffers_.size(); i++)
  {
    fiberPointBuffers_[i].states[0] = -75.0;
    fiberPointBuffers_[i].states[1] = 0.05;
    fiberPointBuffers_[i].states[2] = 0.6;
    fiberPointBuffers_[i].states[3] = 0.325;
  }

}

void FastMonodomainSolver<Control::MultipleInstances<OperatorSplitting::Strang<Control::MultipleInstances<TimeSteppingScheme::Heun<CellmlAdapter<4, 9, FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1> > > > >, Control::MultipleInstances<TimeSteppingScheme::ImplicitEuler<SpatialDiscretization::FiniteElementMethod<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1>, Quadrature::Gauss<2>, Equation::Dynamic::IsotropicDiffusion> > > > > >::
run()
{
  LOG(FATAL) << "not implemented";
  //nestedSolvers_.run();
}

void FastMonodomainSolver<Control::MultipleInstances<OperatorSplitting::Strang<Control::MultipleInstances<TimeSteppingScheme::Heun<CellmlAdapter<4, 9, FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1> > > > >, Control::MultipleInstances<TimeSteppingScheme::ImplicitEuler<SpatialDiscretization::FiniteElementMethod<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1>, Quadrature::Gauss<2>, Equation::Dynamic::IsotropicDiffusion> > > > > >::
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
  std::vector<NestedSolversType::TimeSteppingSchemeType> &instances = nestedSolvers_.instancesLocal();

  for (int i = 0; i < instances.size(); i++)
  {
    instances[i].timeStepping2().writeOutput(0, currentTime_);
  }
}

void FastMonodomainSolver<Control::MultipleInstances<OperatorSplitting::Strang<Control::MultipleInstances<TimeSteppingScheme::Heun<CellmlAdapter<4, 9, FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1> > > > >, Control::MultipleInstances<TimeSteppingScheme::ImplicitEuler<SpatialDiscretization::FiniteElementMethod<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1>, Quadrature::Gauss<2>, Equation::Dynamic::IsotropicDiffusion> > > > > >::
fetchFiberData()
{
  LOG(TRACE) << "fetchFiberData";
  std::vector<NestedSolversType::TimeSteppingSchemeType> &instances = nestedSolvers_.instancesLocal();

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
        
        elementLengthsReceiveBuffer = fiberData_[fiberDataNo].elementLengths.data();
        vmValuesReceiveBuffer = fiberData_[fiberDataNo].vmValues.data();
      }

      for (int rankNo = 0; rankNo < rankSubset->size(); rankNo++)
      {
        nElementsOnRanks[rankNo] = fiberFunctionSpace->meshPartition()->nNodesLocalWithGhosts(0, rankNo) - 1;
        offsetsOnRanks[rankNo] = fiberFunctionSpace->meshPartition()->beginNodeGlobalNatural(0, rankNo);
        nDofsOnRanks[rankNo] = fiberFunctionSpace->meshPartition()->nNodesLocalWithoutGhosts(0, rankNo);
      }

        /*
int MPI_Gatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                void *recvbuf, const int *recvcounts, const int *displs,
                MPI_Datatype recvtype, int root, MPI_Comm comm)
*/
      LOG(DEBUG) << "Gatherv of element lengths to rank " << computingRank << ", values " << localLengths << ", sizes: " << nElementsOnRanks << ", offsets: " << offsetsOnRanks;

      MPI_Gatherv(localLengths.data(), fiberFunctionSpace->nElementsLocal(), MPI_DOUBLE,
                  elementLengthsReceiveBuffer, nElementsOnRanks.data(), offsetsOnRanks.data(),
                  MPI_DOUBLE, computingRank, mpiCommunicator);

      // communicate Vm values
      std::vector<double> vmValuesLocal;
      innerInstances[j].data().solution()->getValuesWithoutGhosts(0, vmValuesLocal);

      LOG(DEBUG) << "Gatherv of values to rank " << computingRank << ", sizes: " << nDofsOnRanks << ", offsets: " << offsetsOnRanks << ", local values " << vmValuesLocal;

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

void FastMonodomainSolver<Control::MultipleInstances<OperatorSplitting::Strang<Control::MultipleInstances<TimeSteppingScheme::Heun<CellmlAdapter<4, 9, FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1> > > > >, Control::MultipleInstances<TimeSteppingScheme::ImplicitEuler<SpatialDiscretization::FiniteElementMethod<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1>, Quadrature::Gauss<2>, Equation::Dynamic::IsotropicDiffusion> > > > > >::
updateFiberData()
{
  // copy Vm from compute buffers to fiberData_
  for (int fiberDataNo = 0; fiberDataNo < fiberData_.size(); fiberDataNo++)
  {
    int nValues = fiberData_[fiberDataNo].vmValues.size();

    for (int valueNo = 0; valueNo < nValues; valueNo++)
    {
      global_no_t valueIndexAllFibers = fiberData_[fiberDataNo].valuesOffset + valueNo;

      global_no_t pointBuffersNo = valueIndexAllFibers / Vc::double_v::Size;
      int entryNo = valueIndexAllFibers % Vc::double_v::Size;

      fiberData_[fiberDataNo].vmValues[valueNo] = fiberPointBuffers_[pointBuffersNo].states[0][entryNo];
    }
  }

  LOG(TRACE) << "updateFiberData";
  std::vector<NestedSolversType::TimeSteppingSchemeType> &instances = nestedSolvers_.instancesLocal();

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

        /*
int MPI_Scatterv(const void *sendbuf, const int *sendcounts, const int *displs, MPI_Datatype sendtype,
                 void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
*/
      // communicate Vm values
      std::vector<double> vmValuesLocal(fiberFunctionSpace->nDofsLocalWithoutGhosts());
      MPI_Scatterv(sendBufferVmValues, nDofsOnRanks.data(), offsetsOnRanks.data(), MPI_DOUBLE,
                   vmValuesLocal.data(), fiberFunctionSpace->nDofsLocalWithoutGhosts(), MPI_DOUBLE,
                   computingRank, mpiCommunicator);

      LOG(DEBUG) << "Scatterv from rank " << computingRank << ", sizes: " << nDofsOnRanks << ", offsets: " << offsetsOnRanks; << ", received local values " << vmValuesLocal

      LOG(DEBUG) << "fiber " << fiberDataNo << ", set values " << vmValuesLocal;
      innerInstances[j].data().solution()->setValuesWithoutGhosts(0, vmValuesLocal);
      instances[i].timeStepping2().instancesLocal()[j].data().solution()->setValuesWithoutGhosts(0, vmValuesLocal);

      // increase index for fiberData_ struct
      if (computingRank == rankSubset->ownRankNo())
        fiberDataNo++;
    }
  }
}

void FastMonodomainSolver<Control::MultipleInstances<OperatorSplitting::Strang<Control::MultipleInstances<TimeSteppingScheme::Heun<CellmlAdapter<4, 9, FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1> > > > >, Control::MultipleInstances<TimeSteppingScheme::ImplicitEuler<SpatialDiscretization::FiniteElementMethod<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1>, Quadrature::Gauss<2>, Equation::Dynamic::IsotropicDiffusion> > > > > >::
computeMonodomain()
{
  LOG(TRACE) << "computeMonodomain";

  // initialize data vector
  // array of vectorized struct

  // fetch timestep widths and total time span
  std::vector<NestedSolversType::TimeSteppingSchemeType> &instances = nestedSolvers_.instancesLocal();

  TimeSteppingScheme::Heun<CellmlAdapterType> &heun = instances[0].timeStepping1().instancesLocal()[0];
  ImplicitEuler &implicitEuler = instances[0].timeStepping2().instancesLocal()[0];
  double prefactor = implicitEuler.discretizableInTime().data().context().getPythonConfig().getOptionDouble("prefactor", 1.0);

  double startTime = instances[0].startTime();
  double timeStepWidthSplitting = instances[0].timeStepWidth();
  const int nTimeStepsSplitting = instances[0].numberTimeSteps();

  heun.setTimeSpan(startTime, startTime + 0.5 * timeStepWidthSplitting);
  double dt0D = heun.timeStepWidth();
  int nTimeSteps0D = heun.numberTimeSteps();

  implicitEuler.setTimeSpan(startTime, startTime + timeStepWidthSplitting);
  double dt1D = implicitEuler.timeStepWidth();
  int nTimeSteps1D = implicitEuler.numberTimeSteps();

  LOG(DEBUG) << "prefactor: " << prefactor << ", dtSplitting: " << timeStepWidthSplitting << ", n steps: " << nTimeStepsSplitting;
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
  for (int timeStepNo = 0; timeStepNo < nTimeStepsSplitting; timeStepNo++)
  {
    // perform Strang splitting
    double currentTime = startTime + timeStepNo * timeStepWidthSplitting;

    LOG(DEBUG) << "splitting " << timeStepNo << "/" << nTimeStepsSplitting << ", t: " << currentTime;

    // compute midTime once per step to reuse it. [currentTime, midTime=currentTime+0.5*timeStepWidth, currentTime+timeStepWidth]
    double midTime = currentTime + 0.5 * timeStepWidthSplitting;

    // perform splitting
    compute0D(currentTime, dt0D, nTimeSteps0D);
    compute1D(currentTime, dt1D, nTimeSteps1D, prefactor);
    compute0D(midTime,     dt0D, nTimeSteps0D);
  }

  currentTime_ = instances[0].endTime();
}

void FastMonodomainSolver<Control::MultipleInstances<OperatorSplitting::Strang<Control::MultipleInstances<TimeSteppingScheme::Heun<CellmlAdapter<4, 9, FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1> > > > >, Control::MultipleInstances<TimeSteppingScheme::ImplicitEuler<SpatialDiscretization::FiniteElementMethod<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1>, Quadrature::Gauss<2>, Equation::Dynamic::IsotropicDiffusion> > > > > >::
compute0D(double startTime, double timeStepWidth, int nTimeSteps)
{
  LOG(DEBUG) << "compute0D(" << startTime << "), " << nTimeSteps << " time steps";

  using Vc::double_v;

  // constants
  const double constant0 = -75;
  const double constant1 = 1;
  const double constant2 = 0;
  const double constant3 = 120;
  const double constant4 = 36;
  const double constant5 = 0.3;
  const double constant6 = constant0 + 115.000;
  const double constant7 = constant0 - 12.0000;
  const double constant8 = constant0 + 10.6130;

  // Heun scheme:
  // y* = y_n + dt*rhs(y_n)
  // y_n+1 = y_n + 0.5*[rhs(y_n) + rhs(y*)]

  // loop over buffers
  for (global_no_t pointBuffersNo = 0; pointBuffersNo < fiberPointBuffers_.size(); pointBuffersNo++)
  {
    int fiberDataNo = pointBuffersNo * Vc::double_v::Size / fiberData_[0].valuesLength;
    int indexInFiber = pointBuffersNo * Vc::double_v::Size - fiberData_[fiberDataNo].valuesOffset;
    int fiberCenterIndex = fiberData_[fiberDataNo].valuesLength / 2;
    bool currentPointIsInCenter = (abs(indexInFiber - fiberCenterIndex) < Vc::double_v::Size);

    int motorUnitNo = fiberData_[fiberDataNo].motorUnitNo;
    VLOG(3) << "pointBuffersNo: " << pointBuffersNo << ", fiberDataNo: " << fiberDataNo << ", indexInFiber: " << indexInFiber << ", motorUnitNo: " << motorUnitNo;

    // get state values to prevent aliasing inefficiencies for the compiler
    double_v state0 = fiberPointBuffers_[pointBuffersNo].states[0];
    double_v state1 = fiberPointBuffers_[pointBuffersNo].states[1];
    double_v state2 = fiberPointBuffers_[pointBuffersNo].states[2];
    double_v state3 = fiberPointBuffers_[pointBuffersNo].states[3];

    VLOG(3) << "  states [" << state0 << "," << state1 << "," << state2 << "," << state3 << "]";

    // loop over timesteps
    for (int timeStepNo = 0; timeStepNo < nTimeSteps; timeStepNo++)
    {
      // determine if fiber gets stimulated
      double currentTime = startTime + timeStepNo * timeStepWidth;
      int index = round(currentTime * setSpecificStatesCallFrequency_);

      double stimulationStartTime = index*setSpecificStatesCallFrequency_;
      double stimulationDuration = currentTime - stimulationStartTime;

      bool stimulate =
        firingEvents_[index % firingEvents_.size()][motorUnitNo % firingEvents_[index % firingEvents_.size()].size()]
        && stimulationDuration <= setSpecificStatesRepeatAfterFirstCall_
        && currentPointIsInCenter;

      if (stimulate)
      {
        LOG(DEBUG) << "stimulate fiber " << fiberData_[fiberDataNo].fiberNo << ", MU " << motorUnitNo << " at t=" << currentTime;
        LOG(DEBUG) << "  pointBuffersNo: " << pointBuffersNo << ", indexInFiber: " << indexInFiber << ", fiberCenterIndex: " << fiberCenterIndex;
        LOG(DEBUG) << "  motorUnitNo: " << motorUnitNo << " (" << motorUnitNo % firingEvents_[index % firingEvents_.size()].size() << ")";
        LOG(DEBUG) << "  time index: " << index << " (" << index % firingEvents_.size() << ")";
        LOG(DEBUG) << "  stimulationStartTime: " << stimulationStartTime << ", stimulationDuration: " << stimulationDuration;
      }

      // perform one step of the heun scheme

      // compute new rates, rhs(y_n)
      const double_v algebraic1                   = ( - 0.100000*(state0+50.0000))/(exp(- (state0+50.0000)/10.0000) - 1.00000);
      const double_v algebraic5                   =  4.00000*exp(- (state0+75.0000)/18.0000);
      const double_v rate1                        =  algebraic1*(1.00000 - state1) -  algebraic5*state1;
      const double_v algebraic2                   =  0.0700000*exp(- (state0+75.0000)/20.0000);
      const double_v algebraic6                   = 1.00000/(exp(- (state0+45.0000)/10.0000)+1.00000);
      const double_v rate2                        =  algebraic2*(1.00000 - state2) -  algebraic6*state2;
      const double_v algebraic3                   = ( - 0.0100000*(state0+65.0000))/(exp(- (state0+65.0000)/10.0000) - 1.00000);
      const double_v algebraic7                   =  0.125000*exp((state0+75.0000)/80.0000);
      const double_v rate3                        =  algebraic3*(1.00000 - state3) -  algebraic7*state3;
      const double_v algebraic0                   =  constant3*state1*state1*state1*state2*(state0 - constant6);
      const double_v algebraic4                   =  constant4*state3*state3*state3*state3*(state0 - constant7);
      const double_v algebraic8                   =  constant5*(state0 - constant8);
      const double_v rate0                        = - (- constant2+algebraic0+algebraic4+algebraic8)/constant1;

      if (pointBuffersNo == fiberPointBuffers_.size()/2)
        VLOG(2) << "increment: [" << rate0*timeStepWidth << "," << rate1*timeStepWidth << "," << rate2*timeStepWidth << "," << rate3*timeStepWidth << "], dt: " << timeStepWidth;

      // intermediate step
      // compute y* = y_n + dt*rhs(y_n), y_n = state, rhs(y_n) = rate, y* = intermediateState
      double_v intermediateState0 = state0 + timeStepWidth*rate0;
      const double_v intermediateState1 = state1 + timeStepWidth*rate1;
      const double_v intermediateState2 = state2 + timeStepWidth*rate2;
      const double_v intermediateState3 = state3 + timeStepWidth*rate3;

      if (stimulate)
      {
        for (int i = 0; i < std::min(3,(int)Vc::double_v::Size); i++)
        {
          intermediateState0[i] = 20.0;
        }
      }

      if (pointBuffersNo == fiberPointBuffers_.size()/2)
        VLOG(2) << "intermediate solution: [" << intermediateState0 << "," << intermediateState1 << "," << intermediateState2 << "," << intermediateState3 << "]";

      // compute new rates, rhs(y*)
      const double_v intermediateAlgebraic1       = ( - 0.100000*(intermediateState0+50.0000))/(exp(- (intermediateState0+50.0000)/10.0000) - 1.00000);
      const double_v intermediateAlgebraic5       =  4.00000*exp(- (intermediateState0+75.0000)/18.0000);
      const double_v intermediateRate1            =  intermediateAlgebraic1*(1.00000 - intermediateState1) -  intermediateAlgebraic5*intermediateState1;
      const double_v intermediateAlgebraic2       =  0.0700000*exp(- (intermediateState0+75.0000)/20.0000);
      const double_v intermediateAlgebraic6       = 1.00000/(exp(- (intermediateState0+45.0000)/10.0000)+1.00000);
      const double_v intermediateRate2            =  intermediateAlgebraic2*(1.00000 - intermediateState2) -  intermediateAlgebraic6*intermediateState2;
      const double_v intermediateAlgebraic3       = ( - 0.0100000*(intermediateState0+65.0000))/(exp(- (intermediateState0+65.0000)/10.0000) - 1.00000);
      const double_v intermediateAlgebraic7       =  0.125000*exp((intermediateState0+75.0000)/80.0000);
      const double_v intermediateRate3            =  intermediateAlgebraic3*(1.00000 - intermediateState3) -  intermediateAlgebraic7*intermediateState3;
      const double_v intermediateAlgebraic0       =  constant3*intermediateState1*intermediateState1*intermediateState1*intermediateState2*(intermediateState0 - constant6);
      const double_v intermediateAlgebraic4       =  constant4*intermediateState3*intermediateState3*intermediateState3*intermediateState3*(intermediateState0 - constant7);
      const double_v intermediateAlgebraic8       =  constant5*(intermediateState0 - constant8);
      const double_v intermediateRate0            = - (- constant2+intermediateAlgebraic0+intermediateAlgebraic4+intermediateAlgebraic8)/constant1;

      // final step
      // y_n+1 = y_n + 0.5*[rhs(y_n) + rhs(y*)]
      const double_v finalState0 = state0 + 0.5*timeStepWidth*(rate0 + intermediateRate0);
      const double_v finalState1 = state1 + 0.5*timeStepWidth*(rate1 + intermediateRate1);
      const double_v finalState2 = state2 + 0.5*timeStepWidth*(rate2 + intermediateRate2);
      const double_v finalState3 = state3 + 0.5*timeStepWidth*(rate3 + intermediateRate3);

      state0 = finalState0;
      state1 = finalState1;
      state2 = finalState2;
      state3 = finalState3;

      if (stimulate)
      {
        for (int i = 0; i < std::min(3,(int)Vc::double_v::Size); i++)
        {
          state0[i] = 20.0;
        }
      }

      if (pointBuffersNo == fiberPointBuffers_.size()/2)
        VLOG(2) << "resulting solution: [" << state0 << "," << state1 << "," << state2 << "," << state3 << "]";

    }

    // store resulting state values
    fiberPointBuffers_[pointBuffersNo].states[0] = state0;
    fiberPointBuffers_[pointBuffersNo].states[1] = state1;
    fiberPointBuffers_[pointBuffersNo].states[2] = state2;
    fiberPointBuffers_[pointBuffersNo].states[3] = state3;

    VLOG(3) << "-> index " << pointBuffersNo << ", states [" << state0 << "," << state1 << "," << state2 << "," << state3 << "]";
  }

  VLOG(1) << "nFiberPointBuffers: " << fiberPointBuffers_.size();
}

void FastMonodomainSolver<Control::MultipleInstances<OperatorSplitting::Strang<Control::MultipleInstances<TimeSteppingScheme::Heun<CellmlAdapter<4, 9, FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1> > > > >, Control::MultipleInstances<TimeSteppingScheme::ImplicitEuler<SpatialDiscretization::FiniteElementMethod<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1>, Quadrature::Gauss<2>, Equation::Dynamic::IsotropicDiffusion> > > > > >::
compute1D(double startTime, double timeStepWidth, int nTimeSteps, double prefactor)
{
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
}

void FastMonodomainSolver<Control::MultipleInstances<OperatorSplitting::Strang<Control::MultipleInstances<TimeSteppingScheme::Heun<CellmlAdapter<4, 9, FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1> > > > >, Control::MultipleInstances<TimeSteppingScheme::ImplicitEuler<SpatialDiscretization::FiniteElementMethod<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1>, Quadrature::Gauss<2>, Equation::Dynamic::IsotropicDiffusion> > > > > >::
reset()
{
  nestedSolvers_.reset();
}

FastMonodomainSolver<Control::MultipleInstances<OperatorSplitting::Strang<Control::MultipleInstances<TimeSteppingScheme::Heun<CellmlAdapter<4, 9, FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1> > > > >, Control::MultipleInstances<TimeSteppingScheme::ImplicitEuler<SpatialDiscretization::FiniteElementMethod<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1>, Quadrature::Gauss<2>, Equation::Dynamic::IsotropicDiffusion> > > > > >::
Data &FastMonodomainSolver<Control::MultipleInstances<OperatorSplitting::Strang<Control::MultipleInstances<TimeSteppingScheme::Heun<CellmlAdapter<4, 9, FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1> > > > >, Control::MultipleInstances<TimeSteppingScheme::ImplicitEuler<SpatialDiscretization::FiniteElementMethod<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1>, Quadrature::Gauss<2>, Equation::Dynamic::IsotropicDiffusion> > > > > >::
data()
{
  return nestedSolvers_.data();
}

void FastMonodomainSolver<Control::MultipleInstances<OperatorSplitting::Strang<Control::MultipleInstances<TimeSteppingScheme::Heun<CellmlAdapter<4, 9, FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1> > > > >, Control::MultipleInstances<TimeSteppingScheme::ImplicitEuler<SpatialDiscretization::FiniteElementMethod<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1>, Quadrature::Gauss<2>, Equation::Dynamic::IsotropicDiffusion> > > > > >::
setTimeSpan(double startTime, double endTime)
{
  nestedSolvers_.setTimeSpan(startTime, endTime);
}

FastMonodomainSolver<Control::MultipleInstances<OperatorSplitting::Strang<Control::MultipleInstances<TimeSteppingScheme::Heun<CellmlAdapter<4, 9, FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1> > > > >, Control::MultipleInstances<TimeSteppingScheme::ImplicitEuler<SpatialDiscretization::FiniteElementMethod<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1>, Quadrature::Gauss<2>, Equation::Dynamic::IsotropicDiffusion> > > > > >::
OutputConnectorDataType FastMonodomainSolver<Control::MultipleInstances<OperatorSplitting::Strang<Control::MultipleInstances<TimeSteppingScheme::Heun<CellmlAdapter<4, 9, FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1> > > > >, Control::MultipleInstances<TimeSteppingScheme::ImplicitEuler<SpatialDiscretization::FiniteElementMethod<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1>, Quadrature::Gauss<2>, Equation::Dynamic::IsotropicDiffusion> > > > > >::
getOutputConnectorData()
{
  return nestedSolvers_.getOutputConnectorData();
}
