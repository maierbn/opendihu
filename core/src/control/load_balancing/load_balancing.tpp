#include "control/load_balancing/load_balancing.h"

#include "spatial_discretization/finite_element_method/finite_element_method.h"

#include <omp.h>
#include <sstream>

namespace Control
{

template<typename CellMLAdapter, typename DiffusionTimeStepping>
LoadBalancing<OperatorSplitting::Strang<TimeSteppingScheme::HeunAdaptiv<CellMLAdapter>,DiffusionTimeStepping>>::
LoadBalancing(DihuContext context) :
  LoadBalancingBase<OperatorSplitting::Strang<TimeSteppingScheme::HeunAdaptiv<CellMLAdapter>,DiffusionTimeStepping>>(context)
{
    // get python config
    this->specificSettings_ = this->context_.getPythonConfig();

    // Set rebalance counter initial to 0. Triggers rebalancing in first step
    this->rebalanceCounter_ = 0;

    // Set rebalancing frequency. In ms of simulation time
    this->rebalanceFrequency_ = 1;
}

template<typename CellMLAdapter, typename DiffusionTimeStepping>
void LoadBalancing<OperatorSplitting::Strang<TimeSteppingScheme::HeunAdaptiv<CellMLAdapter>,DiffusionTimeStepping>>::
rebalance()
{
  // Keyword to find entries in LOG
  LOG(TRACE) << "Rebalance";

  // MPI barrier, to wait for all ranks to reach this code
  std::shared_ptr<Partition::RankSubset> rankSubsetGlobal = this->context_.partitionManager()->rankSubsetForCollectiveOperations();
  MPI_Barrier(rankSubsetGlobal->mpiCommunicator());

  // Get information about finite element method object
  // First, define types
  typedef typename DiffusionTimeStepping::DiscretizableInTime_Type FiniteElementMethodType;
  typedef typename FiniteElementMethodType::FunctionSpace FiberFunctionSpaceType;
  typedef typename DiffusionTimeStepping::Data::FieldVariableType DiffusionFieldVariableType;

  // Second, get used ODE stepping, diffusion stepping and finite element method
  TimeSteppingScheme::HeunAdaptiv<CellMLAdapter> &timeSteppingHeun = this->timeSteppingScheme_.timeStepping1();
  DiffusionTimeStepping &timeSteppingDiffusion = this->timeSteppingScheme_.timeStepping2();
  FiniteElementMethodType &finiteElementMethod = timeSteppingDiffusion.discretizableInTime();

  // Get global ranks and local ranks of a fibre
  std::shared_ptr<Partition::RankSubset> rankSubsetFiber = finiteElementMethod.functionSpace()->meshPartition()->rankSubset();
  LOG(DEBUG) << "rankSubsetGlobal: " << *rankSubsetGlobal;
  LOG(DEBUG) << "rankSubsetFiber: " << *rankSubsetFiber;

  //##########################################################################################################
  //###############################  CHECK FOR NEEDED REBALANCING  ###########################################
  //##########################################################################################################

  // Check if rebalancing is required. Depends on simulation progress
  if (rebalanceCounter_ >= timeSteppingHeun.currentHeunTime())
  {
    // No rebalancing needed, return
    return;
  } else {
    // Defined frequency passed, rebalance. Raise counter for further rebalancing
    rebalanceCounter_ = rebalanceCounter_ + rebalanceFrequency_;
    LOG(DEBUG) << "Starting rebalancing process";
  }

  //##########################################################################################################
  //#######################################  END OF SECTION ##################################################
  //##########################################################################################################

  // Number of elements and elements in current progress and whole fibre
  int nNodesLocalWithoutGhosts = finiteElementMethod.functionSpace()->meshPartition()->nNodesLocalWithoutGhosts();
  int nNodesGlobal = finiteElementMethod.functionSpace()->meshPartition()->nNodesGlobal();
  int nElementsLocal = finiteElementMethod.functionSpace()->meshPartition()->nElementsLocal();
  int nElementsGlobal = finiteElementMethod.functionSpace()->meshPartition()->nElementsGlobal();

  // Log info to keep track in Log-File
  LOG(DEBUG) << "finiteElementMethod, nNodes: local: " << nNodesLocalWithoutGhosts << ", global: " << nNodesGlobal
    << ", nElements: local: " << nElementsLocal << ", global: " << nElementsGlobal;

  // Get current solution of diffusion equation
  std::shared_ptr<DiffusionFieldVariableType> diffusionSolution = timeSteppingDiffusion.data().solution();

  // Log diffusion solution
  LOG(DEBUG) << "finiteElementMethod field variable: " << *diffusionSolution;

  // Get splitting components
  CellMLAdapter &cellMLAdapter = timeSteppingHeun.discretizableInTime();

  // Get cellML adapter information
  int nInstances, nIntermediates, nParameters;
  cellMLAdapter.getNumbers(nInstances, nIntermediates, nParameters);

  // Log cellML adapter information
  LOG(DEBUG) << "cellMLAdapter has " << nInstances << " instances";

  // Define type for cellML solution
  typedef typename TimeSteppingScheme::HeunAdaptiv<CellMLAdapter>::Data::FieldVariableType CellMLFieldVariableType;

  // Get cellML solution
  std::shared_ptr<CellMLFieldVariableType> cellMLSolution = timeSteppingHeun.data().solution();

  //Log cellML solution
  LOG(DEBUG) << "cellML field variable: " << *cellMLSolution;

  // Save all local values

  // Save geometry field values of old mesh
  std::vector<Vec3> geometryFieldValues;
  finiteElementMethod.functionSpace()->geometryField().getValuesWithoutGhosts(geometryFieldValues);

  // Log saved field values of old mesh
  LOG(DEBUG) << "finiteElementMethod has geometryFieldValues: " << geometryFieldValues;

  // Save old diffusion values
  std::vector<double> diffusionValues;
  diffusionSolution->getValuesWithoutGhosts(diffusionValues);

  // Save old cellML values
  const int nCellMLComponents = CellMLAdapter::nComponents();
  std::array<std::vector<double>,nCellMLComponents> cellmlValues;
  cellMLSolution->getValuesWithoutGhosts(cellmlValues);

  //##########################################################################################################
  //#################################  DETERMINE PEAK POSITIONS  #############################################
  //##########################################################################################################

  // Counter to determine amount of peaks in current process
  int peaks_per_process = 0;

  // Vector of peak positions in current process
  std::vector<int> peak_positions_process;

  // Previous peak position. Needed to avoid finding two close points within one stimulus
  int peak_previous = -5;

  // Iterate over vector to determine peak positions of each process
  for(std::vector<double>::iterator it = cellmlValues[0].begin(); it != cellmlValues[0].end(); ++it) {
    if (static_cast<double>(*it) >= 31.0){
        // Check if this peak has already been detected
        if (abs(distance(cellmlValues[0].begin(), it) - peak_previous) > 5) {
            peaks_per_process++;
            peak_positions_process.push_back (distance(cellmlValues[0].begin(), it));
            LOG(DEBUG) << "Peak found at position " << distance(cellmlValues[0].begin(), it);
            peak_previous = distance(cellmlValues[0].begin(), it);
        }
    }
  }

  // Log result of peak detection
  if (peaks_per_process == 0){
    LOG(DEBUG) << "No stimulus found on this process";
  } else {
    LOG(DEBUG) << peaks_per_process <<" peak(s) found on this process";
  }

  //##########################################################################################################
  //#######################################  END OF SECTION ##################################################
  //##########################################################################################################

  //##########################################################################################################
  //##############################  COMMUNICATE PEAK INFORMATION  ############################################
  //##########################################################################################################

  // Gather amount of peaks per fibre. First process of fibre collects the information
  int sum_peaks_fibre = 0;
  MPI_Reduce(&peaks_per_process, &sum_peaks_fibre, 1, MPI_INT, MPI_SUM, 0, rankSubsetFiber->mpiCommunicator());

  // Log info. Only gathering processes log info
  if (rankSubsetFiber->ownRankNo() == 0){
    LOG(DEBUG) << "Gathered infos of " << sum_peaks_fibre <<" peak(s) found on fibre";
  }

  // Gather amount of peaks of all processes. First process (process 0) collects the information
  int sum_peaks_global = 0;
  MPI_Reduce(&peaks_per_process, &sum_peaks_global, 1, MPI_INT, MPI_SUM, 0, rankSubsetGlobal->mpiCommunicator());

  // Log info. Only gathering process logs info
  if (rankSubsetGlobal->ownRankNo() == 0){
    LOG(DEBUG) << "Gathered infos of " << sum_peaks_global <<" peak(s) found globally";
  }

  // Determine maximal amount of peaks per process per fibre. First process of fibre collects the information
  int max_peaks_fibre = 0;
  MPI_Reduce(&peaks_per_process, &max_peaks_fibre, 1, MPI_INT, MPI_MAX, 0, rankSubsetFiber->mpiCommunicator());

  // Log info. Only gathering processes log info
  if (rankSubsetFiber->ownRankNo() == 0){
    LOG(DEBUG) << "Maximum of " << max_peaks_fibre <<" peak(s) per process to this fibre";
  }

  // Determine maximal amount of peaks per process. Process 0 collects the information
  int max_peaks_global = 0;
  MPI_Reduce(&peaks_per_process, &max_peaks_global, 1, MPI_INT, MPI_MAX, 0, rankSubsetGlobal->mpiCommunicator());

  // Log info. Only gathering process logs info
  if (rankSubsetGlobal->ownRankNo() == 0){
    LOG(DEBUG) << "Maximum of " << max_peaks_global <<" peak(s) per process found globally";
  }

  // Braodcast maximal process peak number to all processes of fibre. Needed to define length of send message
  MPI_Bcast(&max_peaks_fibre, 1, MPI_INT, 0, rankSubsetFiber->mpiCommunicator());

  // Log info. All processes are logging
  LOG(DEBUG) << "Received info, that a maximum of " << max_peaks_fibre << " peak(s) are on this fibre per process";

  // Braodcast maximal peak number per process to all processes globally. Needed to define length of send message
  MPI_Bcast(&max_peaks_global, 1, MPI_INT, 0, rankSubsetGlobal->mpiCommunicator());

  // Log info. All processes are logging
  LOG(DEBUG) << "Received info, that a maximum of " << max_peaks_global << " peak(s) are per process globally";

  //Transform local peak position coordinates in global peak position coordinates. Each process does this
  std::vector<int> peak_positions_process_transformed;
  for(std::vector<int>::iterator it = peak_positions_process.begin(); it != peak_positions_process.end(); ++it) {
    peak_positions_process_transformed.push_back(finiteElementMethod.functionSpace()->meshPartition()->getNodeNoGlobalPetsc(static_cast<int>(*it)));
  }

  // Log info. Each process logs its change vector
  LOG(DEBUG) << "Vector of peak positions in global coordinates " << peak_positions_process_transformed;

  //Adjust vector length to maximum number of peaks per process. Missing entries are filled with -1
  while (peak_positions_process_transformed.size() < max_peaks_fibre){
    peak_positions_process_transformed.push_back(-1);
  }

  // Log info
  LOG(DEBUG) << "Vector of peak positions in global coordinates after adjustment " << peak_positions_process_transformed;

  // Collect all peak positions of fibre. First process of fibre collects the information
  std::vector<int> peak_positions_fibre;
  peak_positions_fibre.resize(rankSubsetFiber->size() * max_peaks_fibre);
  MPI_Gather(peak_positions_process_transformed.data(), max_peaks_fibre, MPI_INT, peak_positions_fibre.data(), max_peaks_global, MPI_INT, 0, rankSubsetFiber->mpiCommunicator());

  //Sort vector with peak positions of fibre
  std::sort(peak_positions_fibre.begin(), peak_positions_fibre.end());

  // Log info. Only gathering processes log info
  if (rankSubsetFiber->ownRankNo() == 0){
    LOG(DEBUG) << "Position of peaks on this fiber are " << peak_positions_fibre;
  }

  //Adjust vector length to maximum number of peaks globally. Missing entries are filled with -1
  while (peak_positions_process_transformed.size() < max_peaks_global){
    peak_positions_process_transformed.push_back(-1);
  }

  // Collect all peak positions globally. First process (process 0) collects the information
  std::vector<int> peak_positions_global;
  peak_positions_global.resize(rankSubsetGlobal->size() * max_peaks_global);
  MPI_Gather(peak_positions_process_transformed.data(), max_peaks_global, MPI_INT, peak_positions_global.data(), max_peaks_global, MPI_INT, 0, rankSubsetGlobal->mpiCommunicator());

  //Sort vector with global peak positions
  std::sort(peak_positions_global.begin(), peak_positions_global.end());

  // Log info. Only gathering process logs info
  if (rankSubsetGlobal->ownRankNo() == 0){
    LOG(DEBUG) << "Global position of peaks are " << peak_positions_global;
  }

  //##########################################################################################################
  //######################################  END OF SECTION  ##################################################
  //##########################################################################################################

  //##########################################################################################################
  //#############################  DETERMINE POSSIBLE SPLIT NODES  ###########################################
  //##########################################################################################################

  // Create vector with new possible parting positions
  std::vector<int> split_positions_fibre;

  // Only collecting process of each fibre computes
  if (rankSubsetFiber->ownRankNo() == 0){

    // Add regularly distributed split points to the vector
    for(int rank = 1;  rank < rankSubsetFiber->size(); rank++) {
        split_positions_fibre.push_back(floor(nNodesGlobal/rankSubsetFiber->size())*rank);
    }

    // Add left and right border of peak to vector. A distance of 5 nodes to teh left and to the rigth are chosen
    for(std::vector<int>::iterator it = peak_positions_fibre.begin(); it != peak_positions_fibre.end(); ++it) {
        if (static_cast<int>(*it) - 10 >= 0){
            split_positions_fibre.push_back(static_cast<int>(*it) - 5);
        }
        if (static_cast<int>(*it) + 10 <= nNodesGlobal){
            split_positions_fibre.push_back(static_cast<int>(*it) + 5);
        }
    }

    // Log info
    LOG(DEBUG) << "Possible split positions of this fiber are " << split_positions_fibre;
  }

  //##########################################################################################################
  //######################################  END OF SECTION  ##################################################
  //##########################################################################################################

  //##########################################################################################################
  //#############################  OPTIMISING SPLIT NODE POSITIONS  ##########################################
  //##########################################################################################################

  // Log info
  if (rankSubsetFiber->ownRankNo() == 0){
    LOG(DEBUG) << "Starting optimization process";
  }

  // Final splitting vector
  std::vector<int> split_positions_fibre_final;

  // Resize it for later broadacasting. Size must equal amount of processes - 1
  split_positions_fibre_final.resize(rankSubsetFiber->size() - 1);

  // Only collecting process of each fibre computes
  if (rankSubsetFiber->ownRankNo() == 0){

    // Vector for the permutation containing 0's and 1's
    std::vector<int> split_positions_fibre_to_test(split_positions_fibre.size(),0);

    // Vector for the permutation containing the split point node numbers
    std::vector<int> split_positions_fibre_modified;

    // Add 1's to vector to create basis for binomial coefficient
    for (int i = split_positions_fibre.size() - (rankSubsetFiber->size() -1); i < split_positions_fibre.size(); i++){
      split_positions_fibre_to_test[i] = 1;
    }

    // Variable to hold maximal interval size with a peak found per permutation
    int max_interval_permutation = 0;

    // Variable to hold minimum of found maximal intervals of all permutations
    int min_interval_fibre = nNodesGlobal;

    // Optimize split positions
    do {
      // Clear old vector
      split_positions_fibre_modified.clear();

      // Transform 1's of vector in node numbers
      for(int y = 0; y < split_positions_fibre_to_test.size(); y++){
        if (split_positions_fibre_to_test[y] == 1){
            split_positions_fibre_modified.push_back(split_positions_fibre[y]);
        }
      }

      // Sort vector
      std::sort(split_positions_fibre_modified.begin(), split_positions_fibre_modified.end());

      for(int i = 0; i < split_positions_fibre_modified.size(); i++){
        // When first point is picked
        if (i == 0){
            for (int j = 0; j < peak_positions_fibre.size(); j++){
              if ((0 < peak_positions_fibre[j]) && (peak_positions_fibre[j] < split_positions_fibre_modified[i])){
                  max_interval_permutation = split_positions_fibre_modified[i];
              }
            }
            // Check if only one point is tested
            if (i == split_positions_fibre_modified.size() - 1){
                for (int k = 0; k < peak_positions_fibre.size(); k++){
                    if ((split_positions_fibre_modified[i] < peak_positions_fibre[k]) && (peak_positions_fibre[k] < nNodesGlobal)){
                        max_interval_permutation = std::max(max_interval_permutation, nNodesGlobal - split_positions_fibre_modified[i]);
                    }
                }
            }
        }
        // When last point is picked
        else if (i == split_positions_fibre_modified.size() - 1){
          for (int j = 0; j < peak_positions_fibre.size(); j++){
            if ((split_positions_fibre_modified[i] < peak_positions_fibre[j]) && (peak_positions_fibre[j] < nNodesGlobal)){
                max_interval_permutation = std::max(max_interval_permutation, nNodesGlobal - split_positions_fibre_modified[i]);
            }
          }
        }
        // When point in between is picked
        else {
          for (int j = 0; j < peak_positions_fibre.size(); j++){
            if ((split_positions_fibre_modified[i] < peak_positions_fibre[j]) && (peak_positions_fibre[j] < split_positions_fibre_modified[i + 1])){
                max_interval_permutation = std::max(max_interval_permutation, split_positions_fibre_modified[i +1] - split_positions_fibre_modified[i]);
            }
          }
        }
      }

      // When better solution has been found
      if (max_interval_permutation < min_interval_fibre){

        // Update best solution
        min_interval_fibre = max_interval_permutation;

        // Safe current configuration
        split_positions_fibre_final.clear();

        for (int x = 0; x < split_positions_fibre_modified.size(); x++){
          split_positions_fibre_final.push_back(split_positions_fibre_modified[x]);
        }
      }

    // Log info
    LOG(DEBUG) << "Tested split position set: " << split_positions_fibre_modified << " . Maximal interval with peak found: " << max_interval_permutation;

    } while (std::next_permutation(split_positions_fibre_to_test.begin(),split_positions_fibre_to_test.end()));

    // Log info
    LOG(DEBUG) << "Best split point set calculated to " << split_positions_fibre_final;
  }

  // Broadcast result to all processes of the fiber
  MPI_Bcast(&split_positions_fibre_final[0], split_positions_fibre_final.size(), MPI_INT, 0, rankSubsetFiber->mpiCommunicator());

  // Log successful broadcast
  if (rankSubsetFiber->ownRankNo() != 0){
    LOG(DEBUG) << "Received calculated best split point set " << split_positions_fibre_final;
  }

  //##########################################################################################################
  //######################################  END OF SECTION  ##################################################
  //##########################################################################################################

  //##########################################################################################################
  //##################################  COMMUNICATE NEW VALUES  ##############################################
  //##########################################################################################################

  //------------------------------------  NEW ELEMENT NUMBER  ------------------------------------------------

  // Helping variable to hold start node index of new partition
  int startNodeOfProcessNew;

  // Helping variable to hold start node index of new partition
  int endNodeOfProcessNew;

  // Helping variable to hold start node index of old partition
  int startNodeOfProcessOld = finiteElementMethod.functionSpace()->meshPartition()->beginNodeGlobalNatural(0,-1);

  // Helping variable to hold start node index of old partition
  int endNodeOfProcessOld = startNodeOfProcessOld + nNodesLocalWithoutGhosts - 1;

  // Amount of new local elements
  int nElementsLocalNew_TEST;

  // If first rank is taken
  if (rankSubsetFiber->ownRankNo() == 0){
      startNodeOfProcessNew = 0;
      endNodeOfProcessNew = split_positions_fibre_final[0] - 1;
  }
  // If last rank is taken
  else if (rankSubsetFiber->ownRankNo() == (rankSubsetFiber->size() - 1)){
      startNodeOfProcessNew = split_positions_fibre_final[split_positions_fibre_final.size() - 1];
      endNodeOfProcessNew = nNodesGlobal - 1;
  }
  // If any other rank is taken
  else {
      startNodeOfProcessNew = split_positions_fibre_final[rankSubsetFiber->ownRankNo() - 1];
      endNodeOfProcessNew = split_positions_fibre_final[rankSubsetFiber->ownRankNo()] - 1;
  }

  // Determine new element number
  if (rankSubsetFiber->ownRankNo() == rankSubsetFiber->size() - 1){
      nElementsLocalNew_TEST = endNodeOfProcessNew - startNodeOfProcessNew;
  } else {
      nElementsLocalNew_TEST = endNodeOfProcessNew - startNodeOfProcessNew + 1;
  }

  // Log info
  LOG(DEBUG) << "Changed element number from " << nElementsLocal << " to " << nElementsLocalNew_TEST;

  // Log info
  LOG(DEBUG) << "Old nodes from " << startNodeOfProcessOld << " to " << endNodeOfProcessOld << " changed to " << startNodeOfProcessNew << " to " << endNodeOfProcessNew;

  //--------------------------------------  VALUES OF NODES  -------------------------------------------------

  // Vector of diffusion values to be send and received
  std::vector<double> diffusionValueToSend;
  std::vector<double> diffusionValueToReceive;
  std::vector<double> diffusionValuesNew_TEST = diffusionValues;

  // Vector of nodes to be send
  std::vector<Vec3> nodeValueToSend;
  std::vector<Vec3> nodeValueToReceive;
  std::vector<Vec3> geometryFieldValuesNew_TEST = geometryFieldValues;
  Vec3 temporalNodeValues;

  // Vector of cellml values to be send
  std::array<std::vector<double>,nCellMLComponents> cellmlValueToSend;
  std::array<std::vector<double>,nCellMLComponents> cellmlValueToReceive;
  std::array<std::vector<double>,nCellMLComponents> cellmlValuesNew_TEST = cellmlValues;

  // Counter for each process to determine send/received elements at the left/right end
  int sendReceiveCountLeft = 0;
  int sendReceiveCountRight = 0;

  // Sending, when starting node of process moved to the right
  if (startNodeOfProcessNew > startNodeOfProcessOld){

    // Counter for the while loop
    int loopCounter = 0;

    // Fill sending vector
    while (loopCounter < (startNodeOfProcessNew - startNodeOfProcessOld)){

        // Diffusion values
        diffusionValueToSend.push_back(diffusionValues[loopCounter]);

        // Geometry field values
        nodeValueToSend.push_back(geometryFieldValues[loopCounter]);

        // cellml values
        for (int i = 0; i < nCellMLComponents; i++){
            cellmlValueToSend[i].push_back(cellmlValues[i][loopCounter]);
        }

        // Update counter for to left side of the process
        sendReceiveCountLeft++;

        // Update counter
        loopCounter++;
    }

    // Log info
    LOG(DEBUG) << "Diffusion values to be send " << diffusionValueToSend;
    LOG(DEBUG) << "Nodes to be send " << nodeValueToSend;
    LOG(DEBUG) << "Cellml values to be send " << cellmlValueToSend;

    // Log info
    LOG(DEBUG) << "Creating new vectors ";

    // Create new vectors without elements to be send
    for (int i = 0; i < sendReceiveCountLeft; i++){

        // Remove diffusion values to be send
        diffusionValuesNew_TEST.erase(diffusionValuesNew_TEST.begin());

        // Remove geometry field values
        geometryFieldValuesNew_TEST.erase(geometryFieldValuesNew_TEST.begin());

        // Remove cellml values
        for (int i = 0; i < nCellMLComponents; i++){
            cellmlValuesNew_TEST[i].erase(cellmlValuesNew_TEST[i].begin());
        }
    }

    // Log info
    LOG(DEBUG) << "New vectors for next splitting created";
    LOG(DEBUG) << "New diffusion vector " << diffusionValuesNew_TEST;
    LOG(DEBUG) << "New node vector" << geometryFieldValuesNew_TEST;
    LOG(DEBUG) << "New cellml vector" << cellmlValuesNew_TEST;

    // Send size of vector for diffusion values to previous process
    MPI_Send(&sendReceiveCountLeft, 1, MPI_INT, rankSubsetFiber->ownRankNo() - 1, 0, rankSubsetFiber->mpiCommunicator());

    // Log info
    LOG(DEBUG) << "Length of vectors send to process " << rankSubsetFiber->ownRankNo() - 1;

    // Send diffusion values to previous process
    MPI_Send(diffusionValueToSend.data(), sendReceiveCountLeft, MPI_DOUBLE, rankSubsetFiber->ownRankNo() - 1, 0, rankSubsetFiber->mpiCommunicator());

    // Send nodes to previous process
    for(int i = 0; i < sendReceiveCountLeft; i++){
        MPI_Send(nodeValueToSend[i].data(), 3, MPI_DOUBLE, rankSubsetFiber->ownRankNo() - 1, 0, rankSubsetFiber->mpiCommunicator());
    }

    // Send cellml values to previous process
    for(int i = 0; i < nCellMLComponents; i++){
        MPI_Send(cellmlValueToSend[i].data(), sendReceiveCountLeft, MPI_DOUBLE, rankSubsetFiber->ownRankNo() - 1, 0, rankSubsetFiber->mpiCommunicator());
    }
  }

  // Receive, when end node of process moved to the right
  if (endNodeOfProcessNew > endNodeOfProcessOld){

    // Receive size of vector for diffusion values from successive process
    MPI_Recv(&sendReceiveCountRight, 1, MPI_INT, rankSubsetFiber->ownRankNo() + 1, 0, rankSubsetFiber->mpiCommunicator(), MPI_STATUS_IGNORE);
    LOG(DEBUG) << "Length vectors received from process " << rankSubsetFiber->ownRankNo() + 1;

    // Resize vector
    diffusionValueToReceive.resize(sendReceiveCountRight);
    for(int i = 0; i < nCellMLComponents; i++){
        cellmlValueToReceive[i].resize(sendReceiveCountRight);
    }

    // Receive diffusion values from successive process
    MPI_Recv(diffusionValueToReceive.data(), sendReceiveCountRight, MPI_DOUBLE, rankSubsetFiber->ownRankNo() + 1, 0, rankSubsetFiber->mpiCommunicator(), MPI_STATUS_IGNORE);

    // Log info
    LOG(DEBUG) << "Received diffusion vector " << diffusionValueToReceive;

    for(int i = 0; i < sendReceiveCountRight; i++){
        MPI_Recv(temporalNodeValues.data(), 3, MPI_DOUBLE, rankSubsetFiber->ownRankNo() + 1, 0, rankSubsetFiber->mpiCommunicator(), MPI_STATUS_IGNORE);
        nodeValueToReceive.push_back(temporalNodeValues);
    }

    // Log info
    LOG(DEBUG) << "Received node vector " << nodeValueToReceive;

    for(int i = 0; i < nCellMLComponents; i++){
        MPI_Recv(&cellmlValueToReceive[i][0], sendReceiveCountRight, MPI_DOUBLE, rankSubsetFiber->ownRankNo() + 1, 0, rankSubsetFiber->mpiCommunicator(), MPI_STATUS_IGNORE);
    }

    // Log info
    LOG(DEBUG) << "Received cellml vector " << cellmlValueToReceive;

    // Build new vectors based on received data
    for (int x = 0; x < sendReceiveCountRight; x++){

        // Diffusion vector
        diffusionValuesNew_TEST.push_back(diffusionValueToReceive[x]);

        // Node vector
        geometryFieldValuesNew_TEST.push_back(nodeValueToReceive[x]);

        // Cellml vector
        for (int i = 0; i < nCellMLComponents; i++){
            cellmlValuesNew_TEST[i].push_back(cellmlValueToReceive[i][x]);
        }
    }
  }


  //---------------------------------------------------------------------------------------------------------------------


  // Sending, when ending node of process moved to the left
  if (endNodeOfProcessNew < endNodeOfProcessOld){

    // Counter for the while loop
    int loopCounter = 0;

    // Fill sending vector
    while (loopCounter < (endNodeOfProcessOld - endNodeOfProcessNew)){

        // Diffusion values
        diffusionValueToSend.push_back(diffusionValues[endNodeOfProcessOld - loopCounter]);

        // Geometry fild values
        nodeValueToSend.push_back(geometryFieldValues[loopCounter]);

        // cellml values
        for (int i = 0; i < nCellMLComponents; i++){
            cellmlValueToSend[i].push_back(cellmlValues[i][loopCounter]);
        }

        // Update counter for the right side of the process
        sendReceiveCountRight++;

        // Update variabel
        loopCounter++;
    }

    // Log info
    LOG(DEBUG) << "Diffusion values to be send " << diffusionValueToSend;
    LOG(DEBUG) << "Nodes to be send " << nodeValueToSend;
    LOG(DEBUG) << "Cellml values to be send " << cellmlValueToSend;

    // Remove entries from old vector
    for (int x = 0; x < sendReceiveCountRight; x++){
        // Remove diffusion values
        diffusionValuesNew_TEST.erase(diffusionValuesNew_TEST.end() - 1);

        // Remove geometry field values
        geometryFieldValuesNew_TEST.erase(geometryFieldValuesNew_TEST.end() - 1);

        // Remove cellml values
        for (int i = 0; i < nCellMLComponents; i++){
            cellmlValuesNew_TEST[i].erase(cellmlValuesNew_TEST[i].end() - 1);
        }
    }

    // Log info
    LOG(DEBUG) << "New vectors created";
    LOG(DEBUG) << "New diffusion vector " << diffusionValuesNew_TEST;
    LOG(DEBUG) << "New node vector" << geometryFieldValuesNew_TEST;
    LOG(DEBUG) << "New cellml vector" << cellmlValuesNew_TEST;

    // Send size of vector for diffusion values to previous process
    MPI_Send(&sendReceiveCountRight, 1, MPI_INT, rankSubsetFiber->ownRankNo() + 1, 0, rankSubsetFiber->mpiCommunicator());
    LOG(DEBUG) << "Length of vectors send to process " << rankSubsetFiber->ownRankNo() + 1;

     // Send diffusion values to previous process
    MPI_Send(diffusionValueToSend.data(), sendReceiveCountRight, MPI_DOUBLE, rankSubsetFiber->ownRankNo() + 1, 0, rankSubsetFiber->mpiCommunicator());
    LOG(DEBUG) << "Vector of diffusion values send to process " << rankSubsetFiber->ownRankNo() + 1;

    // Send nodes to previous process
    for(int i = 0; i < sendReceiveCountRight; i++){
        MPI_Send(nodeValueToSend[i].data(), 3, MPI_DOUBLE, rankSubsetFiber->ownRankNo() + 1, 0, rankSubsetFiber->mpiCommunicator());
    }

    // Log info
    LOG(DEBUG) << "Vector of nodes send to process " << rankSubsetFiber->ownRankNo() + 1;

    // Send cellml values to previous process
    for(int i = 0; i < nCellMLComponents; i++){
        MPI_Send(cellmlValueToSend[i].data(), sendReceiveCountRight, MPI_DOUBLE, rankSubsetFiber->ownRankNo() + 1, 0, rankSubsetFiber->mpiCommunicator());
    }

    // Log info
    LOG(DEBUG) << "Vector of cellml values send to process " << rankSubsetFiber->ownRankNo() + 1;
  }

  // Receive, when start node of process moved to the right
  if (startNodeOfProcessNew < startNodeOfProcessOld){

    // Receive size of vector for diffusion values from successive process
    MPI_Recv(&sendReceiveCountLeft, 1, MPI_INT, rankSubsetFiber->ownRankNo() - 1, 0, rankSubsetFiber->mpiCommunicator(), MPI_STATUS_IGNORE);
    LOG(DEBUG) << "Length of vectors received from process " << rankSubsetFiber->ownRankNo() - 1;

    // Resize vectors
    diffusionValueToReceive.resize(sendReceiveCountLeft);
    for(int i = 0; i < nCellMLComponents; i++){
        cellmlValueToReceive[i].resize(sendReceiveCountLeft);
    }

    // Receive diffusion values from successive process
    MPI_Recv(&diffusionValueToReceive[0], sendReceiveCountLeft, MPI_DOUBLE, rankSubsetFiber->ownRankNo() - 1, 0, rankSubsetFiber->mpiCommunicator(), MPI_STATUS_IGNORE);

    // Log info
    LOG(DEBUG) << "Received send diffusion vector " << diffusionValueToReceive;

    // Receive nodes from successive process
    for(int i = 0; i < sendReceiveCountLeft; i++){
        MPI_Recv(temporalNodeValues.data(), 3, MPI_DOUBLE, rankSubsetFiber->ownRankNo() - 1, 0, rankSubsetFiber->mpiCommunicator(), MPI_STATUS_IGNORE);
        nodeValueToReceive.push_back(temporalNodeValues);
    }

    // Log info
    LOG(DEBUG) << "Received send node vector " << nodeValueToReceive;

    // Receive cellml values from successive process
    for(int i = 0; i < nCellMLComponents; i++){
        MPI_Recv(&cellmlValueToReceive[i][0], sendReceiveCountLeft, MPI_DOUBLE, rankSubsetFiber->ownRankNo() - 1, 0, rankSubsetFiber->mpiCommunicator(), MPI_STATUS_IGNORE);
    }

    // Log info
    LOG(DEBUG) << "Received cellml vector " << cellmlValueToReceive;

    // Build new vectors based on received data
    for (int x = 0; x < sendReceiveCountLeft; x++){

        // Define iterators to insert at first psoition
        auto iteratorDiffusion = diffusionValuesNew_TEST.begin();
        auto iteratorNodes = geometryFieldValuesNew_TEST.begin();

        // Diffusion vector
        iteratorDiffusion = diffusionValuesNew_TEST.insert(iteratorDiffusion, diffusionValueToReceive[diffusionValueToReceive.size() -1 - x]);

        // Node vector
        iteratorNodes = geometryFieldValuesNew_TEST.insert(iteratorNodes, nodeValueToReceive[nodeValueToReceive.size() -1 - x]);

        // Cellml vector
        for (int i = 0; i < nCellMLComponents; i++){
            auto iteratorCellml = cellmlValuesNew_TEST[i].begin();
            iteratorCellml = cellmlValuesNew_TEST[i].insert(iteratorCellml, cellmlValueToReceive[i][cellmlValueToReceive[i].size() -1 - x]);
        }
    }
  }

  // Log new infos
  LOG(DEBUG) << " New number of nodes" << endNodeOfProcessNew - startNodeOfProcessNew + 1;
  LOG(DEBUG) << " New element number " << nElementsLocalNew_TEST;
  LOG(DEBUG) << " New diffusion values vector " << diffusionValuesNew_TEST;
  LOG(DEBUG) << " New geometry field values vector " << geometryFieldValuesNew_TEST;
  LOG(DEBUG) << " New cellml value vector " << cellmlValuesNew_TEST;

  LOG(DEBUG) << "YYY" << endNodeOfProcessNew - startNodeOfProcessNew + 1 << ";" << nElementsLocalNew_TEST;

  // Number of elements on own rank after rebalancing
  //int nElementsLocalNew = nElementsLocal;
  int nElementsLocalNew = nElementsLocalNew_TEST;

  // Create new elements based on calculated number of elements
  std::array<int,1> nElementsPerDimensionLocal = {nElementsLocalNew};
  std::array<global_no_t,1> nElementsPerDimensionGlobal;
  std::array<int,1> nRanks = {rankSubsetFiber->size()};

  // Create meshPartition for function space with the currently used ranks
  this->context_.partitionManager()->setRankSubsetForNextCreatedPartitioning(rankSubsetFiber);

  std::shared_ptr<Partition::MeshPartition<FiberFunctionSpaceType>> meshPartition
    = this->context_.partitionManager()->template createPartitioningStructuredLocal<FiberFunctionSpaceType>(
        nElementsPerDimensionGlobal, nElementsPerDimensionLocal, nRanks);

  // Number of nodes on own rank after rebalancing
  int nNodesLocalWithoutGhostsNew = meshPartition->nNodesLocalWithoutGhosts();

  // Vector of node positions
  std::vector<Vec3> nodePositionsWithoutGhosts;
  //nodePositionsWithoutGhosts = geometryFieldValues;
  nodePositionsWithoutGhosts = geometryFieldValuesNew_TEST;
  nodePositionsWithoutGhosts.resize(nNodesLocalWithoutGhostsNew);

  // Used old values again
  //std::vector<double> diffusionValuesNew = diffusionValues;
  std::vector<double> diffusionValuesNew = diffusionValuesNew_TEST;
  //std::array<std::vector<double>,nCellMLComponents> cellmlValuesNew = cellmlValues;
  std::array<std::vector<double>,nCellMLComponents> cellmlValuesNew = cellmlValuesNew_TEST;

  // Log info
  LOG(DEBUG) << "create meshPartition, nElementsPerDimensionGlobal: " << nElementsPerDimensionGlobal;

  // Check for error, when number of global elements changed
  if (nElementsPerDimensionGlobal[0] != nElementsGlobal)
  {
    LOG(FATAL) << "number of global elements changed during rebalancing: old: " << nElementsGlobal << ", new: " << nElementsPerDimensionGlobal[0];
  }

  // Log rebalancing with new values
  LOG(DEBUG) << "rebalancing number of elements: " << nElementsLocal << " -> " << nElementsLocalNew
    << ", number of nodes: " << nNodesLocalWithoutGhosts << " -> " << nNodesLocalWithoutGhostsNew;

  // Create function space with updated mesh partition
  std::stringstream meshName;
  meshName << finiteElementMethod.functionSpace()->meshName() << "_";

  this->context_.partitionManager()->setRankSubsetForNextCreatedPartitioning(rankSubsetFiber);

  std::shared_ptr<FiberFunctionSpaceType> functionSpaceNew = this->context_.meshManager()->template createFunctionSpaceWithGivenMeshPartition<FiberFunctionSpaceType>(
    meshName.str(), meshPartition, nodePositionsWithoutGhosts, nElementsPerDimensionLocal, nRanks);

  finiteElementMethod = std::move(SpatialDiscretization::FiniteElementMethod<
    typename FiberFunctionSpaceType::Mesh,
    typename FiberFunctionSpaceType::BasisFunction,
    typename FiniteElementMethodType::Quadrature,
    typename FiniteElementMethodType::Term
  >(timeSteppingDiffusion.data().context(), functionSpaceNew));


  // save internal value of cellMLAdapter which needs to be preserved
  double lastCallSpecificStatesTime = cellMLAdapter.lastCallSpecificStatesTime();
  LOG(DEBUG) << "store lastCallSpecificStatesTime: " << lastCallSpecificStatesTime;

  // Create new CellMLAdapter with new function space
  cellMLAdapter = std::move(CellMLAdapter(cellMLAdapter, functionSpaceNew));

  // Recreate data structures for cellML
  timeSteppingHeun.reset();
  timeSteppingHeun.initialize();    // retrieves function space from cellMLAdapter

  // set lastCallSpecificStatesTime to previous value
  cellMLAdapter.setLastCallSpecificStatesTime(lastCallSpecificStatesTime);

  // Set all local data
  timeSteppingHeun.data().solution()->setValuesWithoutGhosts(cellmlValuesNew);

  // Recreate data structures for diffusion part
  timeSteppingDiffusion.reset();

  // DIESER BEFEHL MACHT PROBLEME
  timeSteppingDiffusion.initialize();   // retrieves function space from finiteElementMethod

  // Set all local data
  timeSteppingDiffusion.data().solution()->setValuesWithoutGhosts(diffusionValuesNew);

}

} // namespace
