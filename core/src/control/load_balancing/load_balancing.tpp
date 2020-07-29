#include "control/load_balancing/load_balancing.h"

#include "spatial_discretization/finite_element_method/finite_element_method.h"

#include <omp.h>
#include <sstream>

namespace Control
{

template<typename CellMLAdapter, typename DiffusionTimeStepping>
LoadBalancing<OperatorSplitting::Strang<TimeSteppingScheme::HeunAdaptive<CellMLAdapter>,DiffusionTimeStepping>>::
LoadBalancing(DihuContext context) :
  LoadBalancingBase<OperatorSplitting::Strang<TimeSteppingScheme::HeunAdaptive<CellMLAdapter>,DiffusionTimeStepping>>(context)
{
  // get python config
  this->specificSettings_ = this->context_.getPythonConfig();

  // set rebalance counter. 0 triggers rebalancing at beginning of simulation
  this->rebalanceCounter_ = 0;

  // set rebalancing frequency. In ms of simulation time
  this->rebalanceFrequency_ = 1;
}

template<typename CellMLAdapter, typename DiffusionTimeStepping>
void LoadBalancing<OperatorSplitting::Strang<TimeSteppingScheme::HeunAdaptive<CellMLAdapter>,DiffusionTimeStepping>>::
rebalance()
{
  // keyword to find entries in LOG
  LOG(TRACE) << "Rebalance";

  // MPI barrier, to wait for all ranks to reach this code
  std::shared_ptr<Partition::RankSubset> rankSubsetGlobal = this->context_.partitionManager()->rankSubsetForCollectiveOperations();
  MPI_Barrier(rankSubsetGlobal->mpiCommunicator());

  // get information about finite element method object
  // first, define types
  typedef typename DiffusionTimeStepping::DiscretizableInTime FiniteElementMethodType;
  typedef typename FiniteElementMethodType::FunctionSpace FiberFunctionSpaceType;
  typedef typename DiffusionTimeStepping::Data::FieldVariableType DiffusionFieldVariableType;

  // second, get used ODE stepping, diffusion stepping and finite element method
  TimeSteppingScheme::HeunAdaptive<CellMLAdapter> &timeSteppingHeun = this->timeSteppingScheme_.timeStepping1();
  DiffusionTimeStepping &timeSteppingDiffusion = this->timeSteppingScheme_.timeStepping2();
  FiniteElementMethodType &finiteElementMethod = timeSteppingDiffusion.discretizableInTime();

  // get global ranks and local ranks of a fibre
  std::shared_ptr<Partition::RankSubset> rankSubsetFiber = finiteElementMethod.functionSpace()->meshPartition()->rankSubset();
  LOG(DEBUG) << "rankSubsetGlobal: " << *rankSubsetGlobal;
  LOG(DEBUG) << "rankSubsetFiber: " << *rankSubsetFiber;

  // check if rebalancing is required. Depends on simulation progress
  if (this->rebalanceCounter_ >= timeSteppingHeun.currentHeunTime())
  {
    // no rebalancing needed, return
    return;
  }
  else
  {
    // defined frequency passed, rebalance. Raise counter for further rebalancing
    this->rebalanceCounter_ = this->rebalanceCounter_ + rebalanceFrequency_;
    LOG(DEBUG) << "Starting rebalancing process";
  }

  // number of elements and elements in current progress and whole fibre
  int nNodesLocalWithoutGhosts = finiteElementMethod.functionSpace()->meshPartition()->nNodesLocalWithoutGhosts();
  int nNodesGlobal = finiteElementMethod.functionSpace()->meshPartition()->nNodesGlobal();
  int nElementsLocal = finiteElementMethod.functionSpace()->meshPartition()->nElementsLocal();
  int nElementsGlobal = finiteElementMethod.functionSpace()->meshPartition()->nElementsGlobal();

  // log info to keep track in Log-File
  LOG(DEBUG) << "finiteElementMethod, nNodes: local: " << nNodesLocalWithoutGhosts << ", global: " << nNodesGlobal
    << ", nElements: local: " << nElementsLocal << ", global: " << nElementsGlobal;

  // get current solution of diffusion equation
  std::shared_ptr<DiffusionFieldVariableType> diffusionSolution = timeSteppingDiffusion.data().solution();

  // log diffusion solution
  LOG(DEBUG) << "finiteElementMethod field variable: " << *diffusionSolution;

  // get splitting components
  CellMLAdapter &cellMLAdapter = timeSteppingHeun.discretizableInTime();

  // get cellML adapter information
  int nInstances, nAlgebraics, nParameters;
  cellMLAdapter.getNumbers(nInstances, nAlgebraics, nParameters);

  // log cellML adapter information
  LOG(DEBUG) << "cellMLAdapter has " << nInstances << " instances";

  // define type for cellML solution
  typedef typename TimeSteppingScheme::HeunAdaptive<CellMLAdapter>::Data::FieldVariableType CellMLFieldVariableType;

  // get cellML solution
  std::shared_ptr<CellMLFieldVariableType> cellMLSolution = timeSteppingHeun.data().solution();

  // log cellML solution
  LOG(DEBUG) << "cellML field variable: " << *cellMLSolution;

  // save all local values
  // save geometry field values of old mesh
  std::vector<Vec3> geometryFieldValues;
  finiteElementMethod.functionSpace()->geometryField().getValuesWithoutGhosts(geometryFieldValues);

  // log saved field values of old mesh
  LOG(DEBUG) << "finiteElementMethod has geometryFieldValues: " << geometryFieldValues;

  // save old diffusion values
  std::vector<double> diffusionValues;
  diffusionSolution->getValuesWithoutGhosts(diffusionValues);

  // save old cellML values
  const int nCellMLComponents = CellMLAdapter::nStates();
  std::array<std::vector<double>,nCellMLComponents> cellmlValues;
  cellMLSolution->getValuesWithoutGhosts(cellmlValues);

  // counter to determine amount of peaks in current process
  int peaks_per_process = 0;

  // vector of peak positions in current process
  std::vector<int> peak_positions_process;

  // rrevious peak position. Needed to avoid finding two close points within one stimulus
  int peak_previous = -5;

  // iterate over vector to determine peak positions of each process
  for(std::vector<double>::iterator it = cellmlValues[0].begin(); it != cellmlValues[0].end(); ++it)
  {
    if (static_cast<double>(*it) >= 30.0)
    {
      // check if this peak has already been detected
      if (abs(distance(cellmlValues[0].begin(), it) - peak_previous) > 10)
      {
        peaks_per_process++;
        peak_positions_process.push_back (distance(cellmlValues[0].begin(), it));
        LOG(DEBUG) << "Peak found at position " << distance(cellmlValues[0].begin(), it);
        peak_previous = distance(cellmlValues[0].begin(), it);
      }
    }
  }

  // log result of peak detection
  if (peaks_per_process == 0)
  {
    LOG(DEBUG) << "No stimulus found on this process";
  }
  else
  {
    LOG(DEBUG) << peaks_per_process <<" peak(s) found on this process";
  }

  // gather amount of peaks per fibre. First process of fibre collects the information
  int sum_peaks_fibre = 0;
  MPI_Reduce(&peaks_per_process, &sum_peaks_fibre, 1, MPI_INT, MPI_SUM, 0, rankSubsetFiber->mpiCommunicator());
  LOG(INFO) << "There are " << sum_peaks_fibre << " peaks on this fibre";

  // log info. Only gathering processes log info
  if (rankSubsetFiber->ownRankNo() == 0)
  {
    LOG(DEBUG) << "Gathered infos of " << sum_peaks_fibre <<" peak(s) found on fibre";
  }

  // gather amount of peaks of all processes. First process (process 0) collects the information
  int sum_peaks_global = 0;
  MPI_Reduce(&peaks_per_process, &sum_peaks_global, 1, MPI_INT, MPI_SUM, 0, rankSubsetGlobal->mpiCommunicator());

  // log info. Only gathering process logs info
  if (rankSubsetGlobal->ownRankNo() == 0)
  {
    LOG(DEBUG) << "Gathered infos of " << sum_peaks_global <<" peak(s) found globally";
  }

  // determine maximal amount of peaks per process per fibre. First process of fibre collects the information
  int max_peaks_fibre = 0;
  MPI_Reduce(&peaks_per_process, &max_peaks_fibre, 1, MPI_INT, MPI_MAX, 0, rankSubsetFiber->mpiCommunicator());

  // log info. Only gathering processes log info
  if (rankSubsetFiber->ownRankNo() == 0)
  {
    LOG(DEBUG) << "Maximum of " << max_peaks_fibre <<" peak(s) per process to this fibre";
  }

  // determine maximal amount of peaks per process. Process 0 collects the information
  int max_peaks_global = 0;
  MPI_Reduce(&peaks_per_process, &max_peaks_global, 1, MPI_INT, MPI_MAX, 0, rankSubsetGlobal->mpiCommunicator());

  // log info. Only gathering process logs info
  if (rankSubsetGlobal->ownRankNo() == 0)
  {
    LOG(DEBUG) << "Maximum of " << max_peaks_global <<" peak(s) per process found globally";
  }

  // braodcast maximal process peak number to all processes of fibre. Needed to define length of send message
  MPI_Bcast(&max_peaks_fibre, 1, MPI_INT, 0, rankSubsetFiber->mpiCommunicator());

  // log info. All processes are logging
  LOG(DEBUG) << "Received info, that a maximum of " << max_peaks_fibre << " peak(s) are on this fibre per process";

  // braodcast maximal peak number per process to all processes globally. Needed to define length of send message
  MPI_Bcast(&max_peaks_global, 1, MPI_INT, 0, rankSubsetGlobal->mpiCommunicator());

  // log info. All processes are logging
  LOG(DEBUG) << "Received info, that a maximum of " << max_peaks_global << " peak(s) are per process globally";

  // transform local peak position coordinates in global peak position coordinates. Each process does this
  std::vector<int> peak_positions_process_transformed;
  for(std::vector<int>::iterator it = peak_positions_process.begin(); it != peak_positions_process.end(); ++it)
  {
    peak_positions_process_transformed.push_back(finiteElementMethod.functionSpace()->meshPartition()->getNodeNoGlobalPetsc(static_cast<int>(*it)));
  }

  // log info. Each process logs its change vector
  LOG(DEBUG) << "Vector of peak positions in global coordinates " << peak_positions_process_transformed;

  // adjust vector length to maximum number of peaks per process. Missing entries are filled with -1
  while (peak_positions_process_transformed.size() < max_peaks_fibre)
  {
    peak_positions_process_transformed.push_back(-1);
  }

  // log info
  LOG(DEBUG) << "Vector of peak positions in global coordinates after adjustment " << peak_positions_process_transformed;

  // collect all peak positions of fibre. First process of fibre collects the information
  std::vector<int> peak_positions_fibre;
  peak_positions_fibre.resize(rankSubsetFiber->size() * max_peaks_fibre);
  MPI_Gather(peak_positions_process_transformed.data(), max_peaks_fibre, MPI_INT, peak_positions_fibre.data(), max_peaks_global, MPI_INT, 0, rankSubsetFiber->mpiCommunicator());

  // sort vector with peak positions of fibre
  std::sort(peak_positions_fibre.begin(), peak_positions_fibre.end());

  // remove -1 entries
  peak_positions_fibre.erase(std::remove(peak_positions_fibre.begin(), peak_positions_fibre.end(), -1), peak_positions_fibre.end());

  // log info. Only gathering processes log info
  if (rankSubsetFiber->ownRankNo() == 0)
  {
    LOG(DEBUG) << "Position of peaks on this fiber are " << peak_positions_fibre;
    LOG(INFO) << "At time " << timeSteppingHeun.currentHeunTime() << ", the position of peaks on this fiber are " << peak_positions_fibre;
  }

  // adjust vector length to maximum number of peaks globally. Missing entries are filled with -1
  while (peak_positions_process_transformed.size() < max_peaks_global)
  {
    peak_positions_process_transformed.push_back(-1);
  }

  // collect all peak positions globally. First process (process 0) collects the information
  std::vector<int> peak_positions_global;
  peak_positions_global.resize(rankSubsetGlobal->size() * max_peaks_global);
  MPI_Gather(peak_positions_process_transformed.data(), max_peaks_global, MPI_INT, peak_positions_global.data(), max_peaks_global, MPI_INT, 0, rankSubsetGlobal->mpiCommunicator());

  // sort vector with global peak positions
  std::sort(peak_positions_global.begin(), peak_positions_global.end());

  // log info. Only gathering process logs info
  if (rankSubsetGlobal->ownRankNo() == 0)
  {
    LOG(DEBUG) << "Global position of peaks are " << peak_positions_global;
  }

  // create vector with new possible parting positions
  std::vector<int> split_positions_fibre;

  // only collecting process of each fibre computes
  if (rankSubsetFiber->ownRankNo() == 0)
  {

    // add regularly distributed split points to the vector
    for(int rank = 1;  rank < rankSubsetFiber->size(); rank++)
    {
      split_positions_fibre.push_back(floor((nNodesGlobal - 1)/rankSubsetFiber->size()) * rank);
    }

    // add left and right border of peak to vector. A distance of 5 nodes to teh left and to the rigth are chosen
    for(std::vector<int>::iterator it = peak_positions_fibre.begin(); it != peak_positions_fibre.end(); ++it)
    {
      if (static_cast<int>(*it) < (nNodesGlobal/2))
      {
        if (static_cast<int>(*it) - ((25/this->rebalanceFrequency_) - 5) > 0)
        {
          split_positions_fibre.push_back(static_cast<int>(*it) - (25/this->rebalanceFrequency_) - 5);
        }
        if (static_cast<int>(*it) + 5 < nNodesGlobal)
        {
          split_positions_fibre.push_back(static_cast<int>(*it) +5);
        }
      }
      else
      {
        if (static_cast<int>(*it) - 5 > 0)
        {
          split_positions_fibre.push_back(static_cast<int>(*it) - 5);
        }
        if (static_cast<int>(*it) + ((25/this->rebalanceFrequency_) + 8) < nNodesGlobal)
        {
          split_positions_fibre.push_back(static_cast<int>(*it) + (25/this->rebalanceFrequency_) + 8);
        }
      }
    }

    // log info
    LOG(DEBUG) << "Possible split positions of this fiber are " << split_positions_fibre;
  }

  // log info
  if (rankSubsetFiber->ownRankNo() == 0)
  {
    LOG(DEBUG) << "Starting optimization process";
  }

  // final splitting positions
  std::vector<int> split_positions_fibre_final;

  // resize it for later broadacasting. Size must equal amount of processes - 1
  split_positions_fibre_final.resize(rankSubsetFiber->size() - 1);

  // only collecting process of each fibre computes
  if (rankSubsetFiber->ownRankNo() == 0)
  {
    // vector for the permutation containing 0's and 1's
    std::vector<int> split_positions_fibre_to_test(split_positions_fibre.size(),0);

    // vector for the permutation containing the split point node numbers
    std::vector<int> split_positions_fibre_modified;

    // add 1's to vector to create basis for binomial coefficient
    for (int i = split_positions_fibre.size() - (rankSubsetFiber->size() -1); i < split_positions_fibre.size(); i++)
    {
      split_positions_fibre_to_test[i] = 1;
    }

    // variable to hold maximal interval size with a peak found per permutation
    int max_interval_permutation = 0;

    // variable to hold minimum of found maximal intervals of all permutations
    int min_interval_fibre = nNodesGlobal;

    // variable to hold maximal interval size without a peak found per permutation
    int max_interval_peakless_permutation = 0;

    // variable to hold minimum of found maximal intervals without peaks of all permutations
    int min_interval_peakless_fibre = nNodesGlobal;

    // optimize split positions
    do {
      // clear old vector
      split_positions_fibre_modified.clear();

      max_interval_permutation = 0;

      max_interval_peakless_permutation = 0;

      // transform 1's of vector in node numbers
      for(int y = 0; y < split_positions_fibre_to_test.size(); y++)
      {
        if (split_positions_fibre_to_test[y] == 1)
        {
          split_positions_fibre_modified.push_back(split_positions_fibre[y]);
        }
      }

      // sort vector
      std::sort(split_positions_fibre_modified.begin(), split_positions_fibre_modified.end());

      // iterate over all possible split points in current permutation
      for(int i = 0; i < split_positions_fibre_modified.size(); i++)
      {
        // when first point is picked
        if (i == 0)
        {
          // iterate over all peak positions
          for (int j = 0; j < peak_positions_fibre.size(); j++)
          {
            // if there is a peak on this interval
            if ((0 < peak_positions_fibre[j]) && (peak_positions_fibre[j] < split_positions_fibre_modified[i]))
            {
              if (split_positions_fibre_modified[i] > max_interval_permutation)
              {
                max_interval_permutation = split_positions_fibre_modified[i];
              }
            }
            else
            {
              max_interval_peakless_permutation = std::max(max_interval_peakless_permutation, split_positions_fibre_modified[i]);
            }
          }
        }
        // when last point is picked. Check in both directions needed
        else if (i == (split_positions_fibre_modified.size() - 1))
        {
          for (int l = 0; l < peak_positions_fibre.size(); l++)
          {
            if ((split_positions_fibre_modified[i] < peak_positions_fibre[l]) && (peak_positions_fibre[l] < nNodesGlobal - 1))
            {
              max_interval_permutation = std::max(max_interval_permutation, nNodesGlobal - 1 - split_positions_fibre_modified[i]);
            }
            else
            {
              max_interval_peakless_permutation = std::max(max_interval_peakless_permutation, nNodesGlobal - 1 - split_positions_fibre_modified[i]);
            }
            if ((split_positions_fibre_modified[i-1] < peak_positions_fibre[l]) && (peak_positions_fibre[l] < split_positions_fibre_modified[i]))
            {
              max_interval_permutation = std::max(max_interval_permutation, split_positions_fibre_modified[i] - split_positions_fibre_modified[i-1]);
            }
            else
            {
              max_interval_peakless_permutation = std::max(max_interval_peakless_permutation, split_positions_fibre_modified[i] - split_positions_fibre_modified[i-1]);
            }
          }
        }
        // when point in between is picked
        else
        {
          for (int k = 0; k < peak_positions_fibre.size(); k++)
          {
            if ((split_positions_fibre_modified[i-1] < peak_positions_fibre[k]) && (peak_positions_fibre[k] < split_positions_fibre_modified[i]))
            {
              max_interval_permutation = std::max(max_interval_permutation, split_positions_fibre_modified[i] - split_positions_fibre_modified[i-1]);
            }
            else
            {
              max_interval_peakless_permutation = std::max(max_interval_peakless_permutation, split_positions_fibre_modified[i] - split_positions_fibre_modified[i-1]);
            }
          }
        }
      }

      // when better solution has been found
      if (max_interval_permutation < min_interval_fibre)
      {
        // update best solution
        min_interval_fibre = max_interval_permutation;

        // update best peakless interval
        min_interval_peakless_fibre = max_interval_peakless_permutation;

        // safe current configuration
        split_positions_fibre_final.clear();

        for (int x = 0; x < split_positions_fibre_modified.size(); x++)
        {
          split_positions_fibre_final.push_back(split_positions_fibre_modified[x]);
        }
      }

      // when equal solution with same peak distribution is found, compare peakless interval
      if (max_interval_permutation == min_interval_fibre)
      {
        if (max_interval_peakless_permutation < min_interval_peakless_fibre)
        {
          // update best solution
          min_interval_fibre = max_interval_permutation;

          // update best peakless interval
          min_interval_peakless_fibre = max_interval_peakless_permutation;

          // safe current configuration
          split_positions_fibre_final.clear();

          for (int x = 0; x < split_positions_fibre_modified.size(); x++)
          {
            split_positions_fibre_final.push_back(split_positions_fibre_modified[x]);
          }
        }
      }

    // log info
    LOG(DEBUG) << "Tested split position set: " << split_positions_fibre_modified << " . Maximal interval with peak found: " << max_interval_permutation;

    }
    while (std::next_permutation(split_positions_fibre_to_test.begin(),split_positions_fibre_to_test.end()));

    // log info
    LOG(DEBUG) << "Best split point set calculated to " << split_positions_fibre_final;
  }

  // broadcast result to all processes of the fiber
  MPI_Bcast(split_positions_fibre_final.data(), rankSubsetFiber->size() - 1, MPI_INT, 0, rankSubsetFiber->mpiCommunicator());

  // log successful broadcast
  if (rankSubsetFiber->ownRankNo() != 0)
  {
    LOG(DEBUG) << "Received calculated best split point set " << split_positions_fibre_final;
  }

  // helping variable to hold start node index of new partition
  int startNodeOfProcessNew;

  // helping variable to hold start node index of new partition
  int endNodeOfProcessNew;

  // helping variable to hold start node index of old partition
  int startNodeOfProcessOld = finiteElementMethod.functionSpace()->meshPartition()->beginNodeGlobalNatural(0,-1);

  // helping variable to hold start node index of old partition
  int endNodeOfProcessOld = startNodeOfProcessOld + nNodesLocalWithoutGhosts - 1;

  // amount of new local elements
  int nElementsLocalNew_TEST;

  // if first rank is taken
  if (rankSubsetFiber->ownRankNo() == 0)
  {
      startNodeOfProcessNew = 0;
      endNodeOfProcessNew = split_positions_fibre_final[0] - 1;
  }
  // if last rank is taken
  else if (rankSubsetFiber->ownRankNo() == (rankSubsetFiber->size() - 1))
  {
      startNodeOfProcessNew = split_positions_fibre_final[split_positions_fibre_final.size() - 1];
      endNodeOfProcessNew = endNodeOfProcessOld;
  }
  // if any other rank is taken
  else
  {
      startNodeOfProcessNew = split_positions_fibre_final[rankSubsetFiber->ownRankNo() - 1];
      endNodeOfProcessNew = split_positions_fibre_final[rankSubsetFiber->ownRankNo()] - 1;
  }

  // determine new element number
  if (rankSubsetFiber->ownRankNo() == (rankSubsetFiber->size() - 1))
  {
    nElementsLocalNew_TEST = endNodeOfProcessNew - startNodeOfProcessNew;
  }
  else
  {
    nElementsLocalNew_TEST = endNodeOfProcessNew - startNodeOfProcessNew + 1;
  }

  // log info
  LOG(DEBUG) << "Changed element number from " << nElementsLocal << " to " << nElementsLocalNew_TEST;

  // log info
  LOG(DEBUG) << "Old nodes index range from " << startNodeOfProcessOld << " : " << endNodeOfProcessOld << " changed to " << startNodeOfProcessNew << " : " << endNodeOfProcessNew;

  // vector of diffusion values to be send and received
  std::vector<double> diffusionValuesNew_TEST = diffusionValues;

  // vector of nodes to be send
  std::vector<Vec3> geometryFieldValuesNew_TEST = geometryFieldValues;

  // vector of cellml values to be send
  std::array<std::vector<double>,nCellMLComponents> cellmlValuesNew_TEST = cellmlValues;

  // sending, when starting node of process moved to the right
  if (startNodeOfProcessNew > startNodeOfProcessOld)
  {

    // counter for elements to send
    int sendCounterLeft = startNodeOfProcessNew - startNodeOfProcessOld;

    // create data structures
    std::vector<double> diffusionValueToSendLeft;
    std::vector<Vec3> nodeValueToSendLeft;
    std::array<std::vector<double>,nCellMLComponents> cellmlValueToSendLeft;

    // counter for the while loop
    int loopCounter = 0;

    // fill sending vector
    while (loopCounter < sendCounterLeft)
    {
      // diffusion values
      diffusionValueToSendLeft.push_back(diffusionValues[loopCounter]);

      // geometry field values
      nodeValueToSendLeft.push_back(geometryFieldValues[loopCounter]);

      // cellml values
      for (int i = 0; i < nCellMLComponents; i++){
          cellmlValueToSendLeft[i].push_back(cellmlValues[i][loopCounter]);
      }

      // update counter
      loopCounter++;
    }

    // remove diffusion values to be send
    diffusionValuesNew_TEST.erase(diffusionValuesNew_TEST.begin(), diffusionValuesNew_TEST.begin() + sendCounterLeft);

    // remove geometry field values
    geometryFieldValuesNew_TEST.erase(geometryFieldValuesNew_TEST.begin(), geometryFieldValuesNew_TEST.begin() + sendCounterLeft);

    // remove cellml values
    for (int i = 0; i < nCellMLComponents; i++)
    {
      cellmlValuesNew_TEST[i].erase(cellmlValuesNew_TEST[i].begin(), cellmlValuesNew_TEST[i].begin() + sendCounterLeft);
    }

    // log info
    LOG(DEBUG) << "Erased old values, because start point moved to the right";
    LOG(DEBUG) << "New diffusion vector after erasing" << diffusionValuesNew_TEST;
    LOG(DEBUG) << "New node vector after erasing" << geometryFieldValuesNew_TEST;
    LOG(DEBUG) << "New cellml vector after erasing" << cellmlValuesNew_TEST;

    // send diffusion values to previous process
    MPI_Send(diffusionValueToSendLeft.data(), sendCounterLeft, MPI_DOUBLE, rankSubsetFiber->ownRankNo() - 1, 0, rankSubsetFiber->mpiCommunicator());

    // send nodes to previous process
    for(int i = 0; i < sendCounterLeft; i++)
    {
      MPI_Send(nodeValueToSendLeft[i].data(), 3, MPI_DOUBLE, rankSubsetFiber->ownRankNo() - 1, 0, rankSubsetFiber->mpiCommunicator());
    }

    // send cellml values to previous process
    for(int i = 0; i < nCellMLComponents; i++)
    {
      MPI_Send(cellmlValueToSendLeft[i].data(), sendCounterLeft, MPI_DOUBLE, rankSubsetFiber->ownRankNo() - 1, 0, rankSubsetFiber->mpiCommunicator());
    }

    // log info
    LOG(DEBUG) << "Sending vectors to process " << rankSubsetFiber->ownRankNo() - 1 << ", because start point moved to the right ";
    LOG(DEBUG) << "Diffusion values to be send " << diffusionValueToSendLeft;
    LOG(DEBUG) << "Nodes to be send " << nodeValueToSendLeft;
    LOG(DEBUG) << "Cellml values to be send " << cellmlValueToSendLeft;
  }

  // receive, when end node of process moved to the right
  if (endNodeOfProcessNew > endNodeOfProcessOld)
  {
    // counter for received data
    int receiveCounterRight = endNodeOfProcessNew - endNodeOfProcessOld;

    // create data structures
    std::vector<double> diffusionValueToReceiveRight;
    std::vector<Vec3> nodeValueToReceiveRight;
    std::array<std::vector<double>,nCellMLComponents> cellmlValueToReceiveRight;
    Vec3 temporalNodeValuesRight;

    // resize vector
    diffusionValueToReceiveRight.resize(receiveCounterRight);
    for(int i = 0; i < nCellMLComponents; i++)
    {
      cellmlValueToReceiveRight[i].resize(receiveCounterRight);
    }

    // receive diffusion values from successive process
    MPI_Recv(diffusionValueToReceiveRight.data(), receiveCounterRight, MPI_DOUBLE, rankSubsetFiber->ownRankNo() + 1, 0, rankSubsetFiber->mpiCommunicator(), MPI_STATUS_IGNORE);

    for(int i = 0; i < receiveCounterRight; i++)
    {
      MPI_Recv(temporalNodeValuesRight.data(), 3, MPI_DOUBLE, rankSubsetFiber->ownRankNo() + 1, 0, rankSubsetFiber->mpiCommunicator(), MPI_STATUS_IGNORE);
      nodeValueToReceiveRight.push_back(temporalNodeValuesRight);
    }

    for(int i = 0; i < nCellMLComponents; i++)
    {
      MPI_Recv(cellmlValueToReceiveRight[i].data(), receiveCounterRight, MPI_DOUBLE, rankSubsetFiber->ownRankNo() + 1, 0, rankSubsetFiber->mpiCommunicator(), MPI_STATUS_IGNORE);
    }

    // log info
    LOG(DEBUG) << "Received vectors from process " << rankSubsetFiber->ownRankNo() + 1;
    LOG(DEBUG) << "Received diffusion vector " << diffusionValueToReceiveRight;
    LOG(DEBUG) << "Received node vector " << nodeValueToReceiveRight;
    LOG(DEBUG) << "Received cellml vector " << cellmlValueToReceiveRight;

    // build new vectors based on received data
    for (int x = 0; x < receiveCounterRight; x++)
    {
      // diffusion vector
      diffusionValuesNew_TEST.push_back(diffusionValueToReceiveRight[x]);

      // node vector
      geometryFieldValuesNew_TEST.push_back(nodeValueToReceiveRight[x]);

      // cellml vector
      for (int i = 0; i < nCellMLComponents; i++)
      {
        cellmlValuesNew_TEST[i].push_back(cellmlValueToReceiveRight[i][x]);
      }
    }

    // log info
    LOG(DEBUG) << "Created new vector based on received data from process " << rankSubsetFiber->ownRankNo() + 1;
    LOG(DEBUG) << "New diffusion vector " << diffusionValuesNew_TEST;
    LOG(DEBUG) << "New geometry field vector " << geometryFieldValuesNew_TEST;
    LOG(DEBUG) << "New Cellml values vector " << cellmlValuesNew_TEST;
  }

  // sending, when ending node of process moved to the left
  if (endNodeOfProcessNew < endNodeOfProcessOld)
  {
    // counter for elements to send
    int sendCounterRight = endNodeOfProcessOld - endNodeOfProcessNew;

    // create data structures
    std::vector<double> diffusionValueToSendRight;
    std::vector<Vec3> nodeValueToSendRight;
    std::array<std::vector<double>,nCellMLComponents> cellmlValueToSendRight;

    // counter for the while loop
    int loopCounter = 0;

    // fill sending vector
    while (loopCounter < sendCounterRight)
    {
      // diffusion values
      diffusionValueToSendRight.push_back(diffusionValues[diffusionValues.size() - 1 - loopCounter]);

      // geometry field values
      nodeValueToSendRight.push_back(geometryFieldValues[geometryFieldValues.size() - 1 - loopCounter]);

      // cellml values
      for (int i = 0; i < nCellMLComponents; i++){
          cellmlValueToSendRight[i].push_back(cellmlValues[i][cellmlValues[i].size() - 1 -loopCounter]);
      }

      // update variable
      loopCounter++;
    }

    // remove entries from old vector
    for (int x = 0; x < sendCounterRight; x++)
    {

      // remove diffusion values
      diffusionValuesNew_TEST.pop_back();

      // remove geometry field values
      geometryFieldValuesNew_TEST.pop_back();

      // remove cellml values
      for (int i = 0; i < nCellMLComponents; i++){
          cellmlValuesNew_TEST[i].pop_back();
      }
    }

    // log info
    LOG(DEBUG) << "Erased old values, because end point moved to the left";
    LOG(DEBUG) << "New diffusion vector " << diffusionValuesNew_TEST;
    LOG(DEBUG) << "New node vector" << geometryFieldValuesNew_TEST;
    LOG(DEBUG) << "New cellml vector" << cellmlValuesNew_TEST;

    // send diffusion values to previous process
    MPI_Send(diffusionValueToSendRight.data(), sendCounterRight, MPI_DOUBLE, rankSubsetFiber->ownRankNo() + 1, 0, rankSubsetFiber->mpiCommunicator());

    // send nodes to previous process
    for(int i = 0; i < sendCounterRight; i++)
    {
      MPI_Send(nodeValueToSendRight[i].data(), 3, MPI_DOUBLE, rankSubsetFiber->ownRankNo() + 1, 0, rankSubsetFiber->mpiCommunicator());
    }

    // send cellml values to previous process
    for(int i = 0; i < nCellMLComponents; i++)
    {
      MPI_Send(cellmlValueToSendRight[i].data(), sendCounterRight, MPI_DOUBLE, rankSubsetFiber->ownRankNo() + 1, 0, rankSubsetFiber->mpiCommunicator());
    }

    // log info
    LOG(DEBUG) << "Sending vectors to process " << rankSubsetFiber->ownRankNo() + 1 << ", because end point moved to the left ";
    LOG(DEBUG) << "Diffusion values to be send " << diffusionValueToSendRight;
    LOG(DEBUG) << "Nodes to be send " << nodeValueToSendRight;
    LOG(DEBUG) << "Cellml values to be send " << cellmlValueToSendRight;

  }

  // receive, when start node of process moved to the left
  if (startNodeOfProcessNew < startNodeOfProcessOld)
  {
    // counter for received data
    int receiveCounterLeft = startNodeOfProcessOld - startNodeOfProcessNew;

    // create data structures
    std::vector<double> diffusionValueToReceiveLeft;
    std::vector<Vec3> nodeValueToReceiveLeft;
    std::array<std::vector<double>,nCellMLComponents> cellmlValueToReceiveLeft;
    Vec3 temporalNodeValuesLeft;

    // resize vectors
    diffusionValueToReceiveLeft.resize(receiveCounterLeft);
    for(int i = 0; i < nCellMLComponents; i++)
    {
      cellmlValueToReceiveLeft[i].resize(receiveCounterLeft);
    }

    // receive diffusion values from successive process
    MPI_Recv(diffusionValueToReceiveLeft.data(), receiveCounterLeft, MPI_DOUBLE, rankSubsetFiber->ownRankNo() - 1, 0, rankSubsetFiber->mpiCommunicator(), MPI_STATUS_IGNORE);

    // receive nodes from successive process
    for(int i = 0; i < receiveCounterLeft; i++)
    {
      MPI_Recv(temporalNodeValuesLeft.data(), 3, MPI_DOUBLE, rankSubsetFiber->ownRankNo() - 1, 0, rankSubsetFiber->mpiCommunicator(), MPI_STATUS_IGNORE);
      nodeValueToReceiveLeft.push_back(temporalNodeValuesLeft);
    }

    // receive cellml values from successive process
    for(int i = 0; i < nCellMLComponents; i++)
    {
      MPI_Recv(cellmlValueToReceiveLeft[i].data(), receiveCounterLeft, MPI_DOUBLE, rankSubsetFiber->ownRankNo() - 1, 0, rankSubsetFiber->mpiCommunicator(), MPI_STATUS_IGNORE);
    }

    // log info
    LOG(DEBUG) << "Received vectors from process " << rankSubsetFiber->ownRankNo() - 1;
    LOG(DEBUG) << "Received diffusion vector " << diffusionValueToReceiveLeft;
    LOG(DEBUG) << "Received node vector " << nodeValueToReceiveLeft;
    LOG(DEBUG) << "Received cellml vector " << cellmlValueToReceiveLeft;

    // build new vectors based on received data
    for (int x = 0; x < receiveCounterLeft; x++)
    {
      // diffusion vector
      diffusionValuesNew_TEST.insert(diffusionValuesNew_TEST.begin(), diffusionValueToReceiveLeft[x]);

      // node vector
      geometryFieldValuesNew_TEST.insert(geometryFieldValuesNew_TEST.begin(), nodeValueToReceiveLeft[x]);

      // cellml vector
      for (int i = 0; i < nCellMLComponents; i++)
      {
        cellmlValuesNew_TEST[i].insert(cellmlValuesNew_TEST[i].begin(), cellmlValueToReceiveLeft[i][x]);
      }
    }

    // log info
    LOG(DEBUG) << "Created new vector based on received data from process " << rankSubsetFiber->ownRankNo() - 1;
    LOG(DEBUG) << "New diffusion vector " << diffusionValuesNew_TEST;
    LOG(DEBUG) << "New geometry field vector " << geometryFieldValuesNew_TEST;
    LOG(DEBUG) << "New Cellml values vector " << cellmlValuesNew_TEST;
  }

  // log info
  LOG(DEBUG) << "New vectors after rebalancing";
  LOG(DEBUG) << "New diffusion vector " << diffusionValuesNew_TEST;
  LOG(DEBUG) << "New geometry field vector " << geometryFieldValuesNew_TEST;
  LOG(DEBUG) << "New Cellml values vector " << cellmlValuesNew_TEST;

  // number of elements on own rank after rebalancing
  int nElementsLocalNew = nElementsLocalNew_TEST;

  // create new elements based on calculated number of elements
  std::array<int,1> nElementsPerDimensionLocal = {nElementsLocalNew};
  std::array<global_no_t,1> nElementsPerDimensionGlobal;
  std::array<int,1> nRanks = {rankSubsetFiber->size()};

  // create meshPartition for function space with the currently used ranks
  this->context_.partitionManager()->setRankSubsetForNextCreatedPartitioning(rankSubsetFiber);

  std::vector<int> rankNos;
  std::shared_ptr<Partition::MeshPartition<FiberFunctionSpaceType>> meshPartition
    = this->context_.partitionManager()->template createPartitioningStructuredLocal<FiberFunctionSpaceType>(
        this->context_.getPythonConfig(), nElementsPerDimensionGlobal, nElementsPerDimensionLocal, nRanks, rankNos);

  // number of nodes on own rank after rebalancing
  int nNodesLocalWithoutGhostsNew = meshPartition->nNodesLocalWithoutGhosts();

  // vector of node positions
  std::vector<Vec3> nodePositionsWithoutGhosts;
  nodePositionsWithoutGhosts = geometryFieldValuesNew_TEST;
  nodePositionsWithoutGhosts.resize(nNodesLocalWithoutGhostsNew);

  // set new values
  std::vector<double> diffusionValuesNew = diffusionValuesNew_TEST;
  std::array<std::vector<double>,nCellMLComponents> cellmlValuesNew = cellmlValuesNew_TEST;

  // log info
  LOG(DEBUG) << "create meshPartition, nElementsPerDimensionGlobal: " << nElementsPerDimensionGlobal;

  // check for error, when number of global elements changed
  if (nElementsPerDimensionGlobal[0] != nElementsGlobal)
  {
    LOG(FATAL) << "number of global elements changed during rebalancing: old: " << nElementsGlobal << ", new: " << nElementsPerDimensionGlobal[0];
  }

  // log rebalancing with new values
  LOG(DEBUG) << "rebalancing number of elements: " << nElementsLocal << " -> " << nElementsLocalNew
    << ", number of nodes: " << nNodesLocalWithoutGhosts << " -> " << nNodesLocalWithoutGhostsNew;

  // create function space with updated mesh partition
  std::stringstream meshName;
  meshName << finiteElementMethod.functionSpace()->meshName() << "_" << this->rebalanceCounter_;

  this->context_.partitionManager()->setRankSubsetForNextCreatedPartitioning(rankSubsetFiber);

  std::shared_ptr<FiberFunctionSpaceType> functionSpaceNew = this->context_.meshManager()->template createFunctionSpaceWithGivenMeshPartition<FiberFunctionSpaceType>(
    meshName.str(), meshPartition, nodePositionsWithoutGhosts, nElementsPerDimensionLocal, nRanks);

  LOG(DEBUG) << "create FiniteElementMethod object";

  finiteElementMethod = std::move(SpatialDiscretization::FiniteElementMethod<
    typename FiberFunctionSpaceType::Mesh,
    typename FiberFunctionSpaceType::BasisFunction,
    typename FiniteElementMethodType::Quadrature,
    typename FiniteElementMethodType::Term
  >(timeSteppingDiffusion.data().context(), functionSpaceNew));

  // save internal value of cellMLAdapter which needs to be preserved
  double lastCallSpecificStatesTime = cellMLAdapter.lastCallSpecificStatesTime();
  LOG(DEBUG) << "store lastCallSpecificStatesTime: " << lastCallSpecificStatesTime;

  // create new CellMLAdapter with new function space
  cellMLAdapter = std::move(CellMLAdapter(cellMLAdapter, functionSpaceNew));

  // recreate data structures for cellML
  timeSteppingHeun.reset();
  timeSteppingHeun.initialize();    // retrieves function space from cellMLAdapter

  // set lastCallSpecificStatesTime to previous value
  cellMLAdapter.setLastCallSpecificStatesTime(lastCallSpecificStatesTime);

  // set all local data
  timeSteppingHeun.data().solution()->setValuesWithoutGhosts(cellmlValuesNew);

  // recreate data structures for diffusion part
  timeSteppingDiffusion.reset();

  // initialize timeStepping
  timeSteppingDiffusion.initialize();   // retrieves function space from finiteElementMethod

  // set all local data
  timeSteppingDiffusion.data().solution()->setValuesWithoutGhosts(diffusionValuesNew);

  LOG(INFO) << "Reached End of Rebalancing";
}

} // namespace
