#include "specialized_solver/multidomain_solver/multidomain_with_fat_solver.h"

#include <Python.h>  // has to be the first included header

#include <iterator>

namespace TimeSteppingScheme
{


template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusionMuscle,typename FiniteElementMethodDiffusionFat>
void MultidomainWithFatSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle,FiniteElementMethodDiffusionFat>::
findSharedNodesBetweenMuscleAndFat()
{
  using MuscleFunctionSpace = typename FiniteElementMethodDiffusionMuscle::FunctionSpace;
  using FatFunctionSpace = typename FiniteElementMethodDiffusionFat::FunctionSpace;
  std::shared_ptr<MuscleFunctionSpace> functionSpaceMuscle = this->finiteElementMethodDiffusion_.functionSpace();
  std::shared_ptr<FatFunctionSpace> functionSpaceFat = finiteElementMethodFat_.functionSpace();

  // determine nodes that are the same on multiple meshes

  // get node positions of the two meshes
  const int nDofsPerNode = MuscleFunctionSpace::nDofsPerNode();
  assert(nDofsPerNode == FatFunctionSpace::nDofsPerNode());

  std::array<std::vector<Vec3>,2> nodePositions;                            // the node positions of the muscle and the fat mesh
  std::array<std::vector<std::pair<Vec3,node_no_t>>,2> nodePositionsNodes;  // the node positions with nodes for every submesh
  std::array<node_no_t,2> nNodesLocalWithoutGhosts = {functionSpaceMuscle->nNodesLocalWithoutGhosts(), functionSpaceFat->nNodesLocalWithoutGhosts()};

  // get geometry fields for function spaces
  functionSpaceMuscle->geometryField().getValuesWithoutGhosts(nodePositions[0]);
  functionSpaceFat->geometryField().getValuesWithoutGhosts(nodePositions[1]);


  // iterate over muscle and fat mesh and save all node positions
  for(int i = 0; i < 2; i++)
  {
    // store node positions with local node nos
    for (node_no_t nodeNoLocal = 0; nodeNoLocal < nNodesLocalWithoutGhosts[i]; nodeNoLocal++)
    {
      nodePositionsNodes[i].push_back(std::make_pair(nodePositions[i][nodeNoLocal*nDofsPerNode + 0], nodeNoLocal));
    }

    // sort according to x coordinate of node positions
    std::sort(nodePositionsNodes[i].begin(), nodePositionsNodes[i].end(), [](const std::pair<Vec3,node_no_t> &a, const std::pair<Vec3,node_no_t> &b)
    {
      return a.first[0] < b.first[0];
    });
  }

  const double nodePositionEqualTolerance = 1e-5;

  // for muscle mesh find shared nodes in fat mesh

  // loop over node positions of muscle mesh
  for (node_no_t nodeNoLocal = 0; nodeNoLocal < nodePositions[0].size(); nodeNoLocal++)
  {
    Vec3 position = nodePositions[0][nodeNoLocal];
    VLOG(2) << "nodeNo " << nodeNoLocal << ", position " << position << " at x=" << position[0];

    // for fat layer mesh
    int indexOtherMesh = 1;
  
    // find node that has closest x coordinate
    VLOG(2) << "  fat mesh, find node that has the closest x coordinate to " << position[0];
    VLOG(2) << "  get last node with x coordinate that is lower than " << position[0] << " by more than tolerance " << nodePositionEqualTolerance;

    // get last node with x coordinate that is lower by more than tolerance
    int k = nodePositionsNodes[indexOtherMesh].size() / 2;
    int kPrevious = -1;
    int lower = 0;
    int upper = nodePositionsNodes[indexOtherMesh].size();

    if (upper > 0)
      while (k != kPrevious)
      {
        Vec3 currentNodePosition = nodePositionsNodes[indexOtherMesh][k].first;
        if (currentNodePosition[0] < position[0]-nodePositionEqualTolerance)
        {
          lower = k;
        }
        else
        {
          upper = k;
        }
        kPrevious = k;
        k = (upper + lower) / 2;

        VLOG(2) << "  range [" << lower << "," << upper << "] k:" << k << ", x:" << currentNodePosition[0];
      }
    VLOG(2) << "  now check all node positions of otherMesh " << indexOtherMesh << " that have x=" << position[0] << " within tolerance";

    // check all node positions after k
    for (;k < nodePositionsNodes[indexOtherMesh].size(); k++)
    {
      Vec3 nodePositionOtherMesh = nodePositionsNodes[indexOtherMesh][k].first;
      node_no_t nodeNoLocalOtherMesh = nodePositionsNodes[indexOtherMesh][k].second;

      assert (nodeNoLocal < functionSpaceMuscle->nNodesLocalWithoutGhosts());
      assert (nodeNoLocalOtherMesh < functionSpaceFat->nNodesLocalWithoutGhosts());

      if (nodePositionOtherMesh[0] > position[0]+nodePositionEqualTolerance)
      {
        VLOG(2) << "  node k: " << k << ", nodeNo: " << nodeNoLocalOtherMesh << ", position: " << nodePositionOtherMesh
          << " is over " << position[0]+nodePositionEqualTolerance << " -> break";
        break;
      }

      double distance = MathUtility::distance<3>(position, nodePositionOtherMesh);
      VLOG(2) << "  node k: " << k << ", nodeNo: " << nodeNoLocalOtherMesh << ", position: " << nodePositionOtherMesh << ", distance: " << distance;

      // if the other mesh node is at the same position as the first node
      if (distance <= nodePositionEqualTolerance)
      {
        VLOG(2) << "   node is shared.";
        sharedNodes_[nodeNoLocal] = nodeNoLocalOtherMesh;

        // add dofs at node to set of border nodes of fat mesh
        for (int i = 0; i < FunctionSpace::nDofsPerNode(); i++)
        {
          dof_no_t dofNoLocalFat = FunctionSpace::nDofsPerNode()*nodeNoLocalOtherMesh + i;
          borderDofsFat_.insert(dofNoLocalFat);
        }

        break;
      }
    }
  }

  nSharedDofsLocal_ = borderDofsFat_.size();

  LOG(DEBUG) << "functionSpaceMuscle: " << *functionSpaceMuscle->meshPartition();
  LOG(DEBUG) << "functionSpaceFat: " << *functionSpaceFat->meshPartition();

  LOG(DEBUG) << sharedNodes_.size() << "sharedNodes_ (\"muscle mesh dof\": fat mesh dof):\n" << sharedNodes_;
  LOG(DEBUG) << borderDofsFat_.size() << " borderDofsFat_:\n" << borderDofsFat_;
  
  LOG(DEBUG) << "n dofs muscle (global): " << functionSpaceMuscle->nDofsGlobal() << " (local without ghosts: " 
    << functionSpaceMuscle->nDofsLocalWithoutGhosts() << ", local with ghosts: " << functionSpaceMuscle->nDofsLocalWithGhosts() << ")";
  LOG(DEBUG) << "n dofs fat (global): " << functionSpaceFat->nDofsGlobal() << " (local without ghosts: " 
    << functionSpaceFat->nDofsLocalWithoutGhosts() << ", local with ghosts: " << functionSpaceFat->nDofsLocalWithGhosts() << ")";
  LOG(DEBUG) << "n shared dofs (global): " << borderDofsGlobalFat_.size() << " (local: " << nSharedDofsLocal_ << ")";


  // globally exchange border dofs
  int nRanks = functionSpaceFat->meshPartition()->rankSubset()->size();
  int ownRankNo = functionSpaceFat->meshPartition()->rankSubset()->ownRankNo();
  MPI_Comm mpiCommunicator =  functionSpaceFat->meshPartition()->rankSubset()->mpiCommunicator();

  std::vector<std::vector<global_no_t>> sharedNodesFatOnRanks(nRanks);   //< for every rank the shared nodes in the fat mesh global petsc numbering
  std::vector<std::vector<PetscInt>> sharedDofsGlobalOnRanks(nRanks);

  // store mapping between fat and muscle border dof nos for own rank
  sharedDofsGlobalOnRanks[ownRankNo].reserve(sharedNodes_.size()*2);    // mapping between fat and muscle shared dofs 

  for (std::pair<node_no_t,node_no_t> sharedNodes : sharedNodes_)
  {
    for (int i = 0; i < FunctionSpace::nDofsPerNode(); i++)
    {
      dof_no_t dofNoLocalFat = sharedNodes.second*FunctionSpace::nDofsPerNode() + i;
      dof_no_t dofNoLocalMuscle = sharedNodes.first*FunctionSpace::nDofsPerNode() + i;

      LOG(DEBUG) << "dofNoLocal fat: " << dofNoLocalFat << ", muscle: " << dofNoLocalMuscle;

      global_no_t dofNoFatGlobalPetsc = functionSpaceFat->meshPartition()->getDofNoGlobalPetsc(dofNoLocalFat);
      global_no_t dofNoMuscleGlobalPetsc = functionSpaceMuscle->meshPartition()->getDofNoGlobalPetsc(dofNoLocalMuscle);

      sharedDofsGlobalOnRanks[ownRankNo].push_back(dofNoFatGlobalPetsc);
      sharedDofsGlobalOnRanks[ownRankNo].push_back(dofNoMuscleGlobalPetsc);
    }
  }
  
  // store global nos from borderDofsFat_ in sharedNodesFatOnRanks[ownRankNo]
  sharedNodesFatOnRanks[ownRankNo].reserve(borderDofsFat_.size());
  for (dof_no_t dofNoLocalFat : borderDofsFat_)
  {
    global_no_t dofNoGlobalPetsc = functionSpaceFat->meshPartition()->getDofNoGlobalPetsc(dofNoLocalFat);
    sharedNodesFatOnRanks[ownRankNo].push_back(dofNoGlobalPetsc);
  }

  // determine number of shared nodes on all the ranks
  std::vector<PetscInt> nSharedNodesOnRanks(nRanks);
  nSharedNodesOnRanks[ownRankNo] = borderDofsFat_.size();

  // communicate how many shared nodes there are on every rank
  MPIUtility::handleReturnValue(MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, nSharedNodesOnRanks.data(),
                                              1, MPIU_INT, mpiCommunicator), "MPI_Allgather");

  // communicate shared dof nos, sharedDofsGlobalOnRanks and shared nodes on fat mesh, sharedNodesFatOnRanks
  PetscInt nEntriesTotal = 0;
  for (int rankNo = 0; rankNo < nRanks; rankNo++)                                            
  {
    int nEntries = nSharedNodesOnRanks[rankNo];
    nEntriesTotal += nEntries;
    sharedDofsGlobalOnRanks[rankNo].resize(nEntries*2);
    sharedNodesFatOnRanks[rankNo].resize(nEntries);

    // broadcast node nos to all ranks
    MPI_Bcast(sharedDofsGlobalOnRanks[rankNo].data(), 2*nEntries, MPIU_INT, rankNo, mpiCommunicator);
    MPI_Bcast(sharedNodesFatOnRanks[rankNo].data(), nEntries, MPI_LONG_LONG_INT, rankNo, mpiCommunicator);
  }

  // add all received global nos to a single set
  for (int rankNo = 0; rankNo < nRanks; rankNo++)                                            
  {
    // add global nos for fat mesh
    borderDofsGlobalFat_.insert(sharedNodesFatOnRanks[rankNo].begin(), sharedNodesFatOnRanks[rankNo].end());

    // add mapping from fat mesh dofs to muscle mesh dofs to a single map
    assert(sharedDofsGlobalOnRanks[rankNo].size() == nSharedNodesOnRanks[rankNo]*2);
    for (int entryNo = 0; entryNo < sharedDofsGlobalOnRanks[rankNo].size(); entryNo += 2)
    {
      fatDofToMuscleDofGlobal_.insert(std::make_pair(sharedDofsGlobalOnRanks[rankNo][entryNo], sharedDofsGlobalOnRanks[rankNo][entryNo+1]));    // (fat mesh dof, muscle mesh dof)
      muscleDofToFatDofGlobal_.insert(std::make_pair(sharedDofsGlobalOnRanks[rankNo][entryNo+1], sharedDofsGlobalOnRanks[rankNo][entryNo]));    // (muscle mesh dof, fat mesh dof)
    }
  }

  // output result
#ifndef NDEBUG
  LOG(DEBUG) << borderDofsGlobalFat_.size() << " borderDofsGlobalFat_:\n" << borderDofsGlobalFat_;
  std::stringstream s;
  for (PetscInt dofNoGlobal = 0; dofNoGlobal < this->finiteElementMethodFat_.data().functionSpace()->nDofsGlobal(); dofNoGlobal++)
  {
    if (borderDofsGlobalFat_.find(dofNoGlobal) != borderDofsGlobalFat_.end())
    {
      s << dofNoGlobal << ":shared, ";
    }
    else
    {
      s << dofNoGlobal << ":" << getDofNoGlobalFatWithoutSharedDofs(dofNoGlobal) << ", ";
    }
  }
  LOG(DEBUG) << "global nos without shared dofs: " << s.str();

  s.str("");
  for (PetscInt dofNoGlobalFat : borderDofsGlobalFat_)
  {
    PetscInt dofNoGlobalMuscle = getDofNoGlobalMuscleFromDofNoGlobalFat(dofNoGlobalFat);
     
    s << dofNoGlobalFat << ":" << dofNoGlobalMuscle << ", ";
  }
  
  LOG(DEBUG) << fatDofToMuscleDofGlobal_.size() << " fatDofToMuscleDofGlobal_: " << fatDofToMuscleDofGlobal_;
  LOG(DEBUG) << muscleDofToFatDofGlobal_.size() << " muscleDofToFatDofGlobal_: " << muscleDofToFatDofGlobal_;
  LOG(DEBUG) << borderDofsGlobalFat_.size() << " borderDofsGlobalFat_: " << borderDofsGlobalFat_;
  LOG(DEBUG) << "border dofs fat:muscle " << s.str();
#endif
}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusionMuscle,typename FiniteElementMethodDiffusionFat>
PetscInt MultidomainWithFatSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle,FiniteElementMethodDiffusionFat>::
getDofNoGlobalFatWithoutSharedDofs(PetscInt dofNoGlobal)
{
  int nSharedBefore = std::distance(borderDofsGlobalFat_.begin(), borderDofsGlobalFat_.upper_bound(dofNoGlobal));
  PetscInt dofNoGlobalWithoutSharedDofs = dofNoGlobal - nSharedBefore;
  return dofNoGlobalWithoutSharedDofs;
}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusionMuscle,typename FiniteElementMethodDiffusionFat>
PetscInt MultidomainWithFatSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle,FiniteElementMethodDiffusionFat>::
getDofNoGlobalMuscleFromDofNoGlobalFat(PetscInt dofNoGlobal)
{
  assert(fatDofToMuscleDofGlobal_.find(dofNoGlobal) != fatDofToMuscleDofGlobal_.end());
  return fatDofToMuscleDofGlobal_[dofNoGlobal];
}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusionMuscle,typename FiniteElementMethodDiffusionFat>
PetscInt MultidomainWithFatSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle,FiniteElementMethodDiffusionFat>::
getDofNoGlobalFatFromDofNoGlobalMuscle(PetscInt dofNoGlobal)
{
  assert(muscleDofToFatDofGlobal_.find(dofNoGlobal) != muscleDofToFatDofGlobal_.end());
  return muscleDofToFatDofGlobal_[dofNoGlobal];
}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusionMuscle,typename FiniteElementMethodDiffusionFat>
void MultidomainWithFatSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle,FiniteElementMethodDiffusionFat>::
initializeBorderVariables()
{
  // system to be solved:
  //
  // [A^1_Vm,Vm   |            |             | B^1_Vm,phie |             ]   [ V^1_m^(i+1) ]    [b^1^(i)]
  // [            | A^2_Vm,Vm  |             | B^2_Vm,phie |             ]   [ V^2_m^(i+1) ]    [b^2^(i)]
  // [   ...      |            | A^M_Vm,Vm   | B^M_Vm,phie |             ] * [ V^M_m^(i+1) ] =  [b^M^(i)]
  // [B^1_phie,Vm |B^2_phie,Vm | B^M_phie,Vm | B_phie,phie |      D      ]   [ phi_e^(i+1) ]    [ 0     ]
  // [            |            |             |     E       | C_phib,phib ]   [ phi_b^(i+1) ]    [ 0     ]

  // This method initializes C, D and E, also the rhs

  PetscErrorCode ierr;
  MPI_Comm mpiCommunicator = this->dataMultidomain_.functionSpace()->meshPartition()->mpiCommunicator();

  // -------
  // create vector for rhs
  // initialize PETSc vector object
  Vec vecZero;
  ierr = VecCreate(mpiCommunicator, &vecZero); CHKERRV(ierr);

  // set name of vector
  std::string name("vecZero");
  ierr = PetscObjectSetName((PetscObject) vecZero, name.c_str()); CHKERRV(ierr);

  // initialize size of vector
  int nEntriesLocal = this->finiteElementMethodFat_.data().functionSpace()->nDofsLocalWithoutGhosts() - nSharedDofsLocal_; 
  ierr = VecSetSizes(vecZero, nEntriesLocal, PETSC_DECIDE); CHKERRV(ierr);

  // set sparsity type and other options
  ierr = VecSetFromOptions(vecZero); CHKERRV(ierr);

  // set all entries to 0.0
  ierr = VecZeroEntries(vecZero); CHKERRV(ierr);

  // store the vector at the bottom of the rhs
  this->subvectorsRightHandSide_[this->nCompartments_+1] = vecZero;

  // -------
  // create vector for solution
  // initialize PETSc vector object
  Vec vecSolutionPhiB;
  ierr = VecCreate(mpiCommunicator, &vecSolutionPhiB); CHKERRV(ierr);

  // set name of vector
  name = "phi_b";
  ierr = PetscObjectSetName((PetscObject) vecSolutionPhiB, name.c_str()); CHKERRV(ierr);

  // initialize size of vector
  nEntriesLocal = this->finiteElementMethodFat_.data().functionSpace()->nDofsLocalWithoutGhosts() - nSharedDofsLocal_; 
  ierr = VecSetSizes(vecSolutionPhiB, nEntriesLocal, PETSC_DECIDE); CHKERRV(ierr);

  // set sparsity type and other options
  ierr = VecSetFromOptions(vecSolutionPhiB); CHKERRV(ierr);

  // set all entries to 0.0
  ierr = VecZeroEntries(vecSolutionPhiB); CHKERRV(ierr);

  // store the vector at the bottom of the solution
  this->subvectorsSolution_[this->nCompartments_+1] = vecSolutionPhiB;
  

  // -------
  // create copy of matrix B_phie,phie
  Mat originalMatrixB = this->finiteElementMethodDiffusionTotal_.data().stiffnessMatrix()->valuesGlobal();
  Mat matrixB;
  //ierr = MatDuplicate(originalMatrixB, MAT_COPY_VALUES, &matrixB); CHKERRV(ierr);
  ierr = MatConvert(originalMatrixB, MATSAME, MAT_INITIAL_MATRIX, &matrixB); CHKERRV(ierr);

  // initialize name
  name = "B_phie,phie";
  ierr = PetscObjectSetName((PetscObject)matrixB, name.c_str()); CHKERRV(ierr);

  // -------
  // create matrix c, which the original C_phib,phib matrix but without rows and columns that correspond to shared dofs between muscle and fat mesh
  Mat originalMatrixC;
  if (enableFatComputation_)
  {
    originalMatrixC = this->finiteElementMethodFat_.data().stiffnessMatrix()->valuesGlobal();
  }
  else
  {
    // if the fat layer should not be computed, this is for debugging if the option "enableFatComputation" is set to False, set originalMatrixC to be the identity
    ierr = MatConvert(this->finiteElementMethodFat_.data().stiffnessMatrix()->valuesGlobal(), MATSAME, MAT_INITIAL_MATRIX, &originalMatrixC); CHKERRV(ierr);
    ierr = MatZeroEntries(originalMatrixC); CHKERRV(ierr);
    ierr = MatShift(originalMatrixC, 1.0); CHKERRV(ierr);
  }

  PetscInt nRows, nColumns;
  ierr = MatGetSize(originalMatrixC, &nRows, &nColumns); CHKERRV(ierr);
  LOG(DEBUG) << "originalMatrixC: " << nRows << "x" << nColumns;
  LOG(DEBUG) << "n dofs muscle (global): " << this->finiteElementMethodDiffusionTotal_.data().functionSpace()->nDofsGlobal() << " (local: " << this->finiteElementMethodDiffusionTotal_.data().functionSpace()->nDofsLocalWithoutGhosts() << ")";
  LOG(DEBUG) << "n dofs fat (global): " << this->finiteElementMethodFat_.data().functionSpace()->nDofsGlobal() << " (local: " << this->finiteElementMethodFat_.data().functionSpace()->nDofsLocalWithoutGhosts() << ")";
  LOG(DEBUG) << "n shared dofs (global): " << borderDofsGlobalFat_.size() << " (local: " << nSharedDofsLocal_ << ")";
  
  // get number of allocated non-zeros of full c matrix
  MatInfo matInfo;
  ierr = MatGetInfo(this->finiteElementMethodFat_.data().stiffnessMatrix()->valuesGlobal(), MAT_GLOBAL_MAX, &matInfo); CHKERRV(ierr);
  PetscInt nNonZerosDiagonal = (PetscInt)matInfo.nz_allocated;
  
  // create matrix
  Mat matrixC;
  ierr = MatCreate(mpiCommunicator, &matrixC); CHKERRV(ierr);

  // initialize size
  int nRowsLocal = this->finiteElementMethodFat_.data().functionSpace()->nDofsLocalWithoutGhosts() - nSharedDofsLocal_;
  int nColumnsLocal = nRowsLocal;
  ierr = MatSetSizes(matrixC, nRowsLocal, nColumnsLocal, PETSC_DECIDE, PETSC_DECIDE); CHKERRV(ierr);
  ierr = MatSetType(matrixC, MATAIJ); CHKERRV(ierr);
  
  // sparse matrix: preallocation of internal data structure
  // http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MATAIJ.html#MATAIJ
  // MATAIJ = "aij" - A matrix type to be used for sparse matrices. This matrix type is identical to MATSEQAIJ when constructed with a single process communicator, and MATMPIAIJ otherwise.
  // As a result, for single process communicators, MatSeqAIJSetPreallocation is supported, and similarly MatMPIAIJSetPreallocation is supported for communicators controlling multiple processes.
  // It is recommended that you call both of the above preallocation routines for simplicity.
  PetscInt nNonZerosOffdiagonal = nNonZerosDiagonal;
  ierr = MatSeqAIJSetPreallocation(matrixC, nNonZerosDiagonal, NULL); CHKERRV(ierr);
  ierr = MatMPIAIJSetPreallocation(matrixC, nNonZerosDiagonal, NULL, nNonZerosOffdiagonal, NULL); CHKERRV(ierr);
  
  // initialize name
  name = "C_phib,phib";
  ierr = PetscObjectSetName((PetscObject)matrixC, name.c_str()); CHKERRV(ierr);

  // --------
  // create matrix D
  Mat matrixD;
  ierr = MatCreate(mpiCommunicator, &matrixD); CHKERRV(ierr);

  // initialize size
  nRowsLocal = this->finiteElementMethodDiffusionTotal_.data().functionSpace()->nDofsLocalWithoutGhosts();
  nColumnsLocal = this->finiteElementMethodFat_.data().functionSpace()->nDofsLocalWithoutGhosts() - nSharedDofsLocal_;
  ierr = MatSetSizes(matrixD, nRowsLocal, nColumnsLocal, PETSC_DECIDE, PETSC_DECIDE); CHKERRV(ierr);
  ierr = MatSetType(matrixD, MATAIJ); CHKERRV(ierr);
  
  ierr = MatSeqAIJSetPreallocation(matrixD, nNonZerosDiagonal, NULL); CHKERRV(ierr);
  ierr = MatMPIAIJSetPreallocation(matrixD, nNonZerosDiagonal, NULL, nNonZerosOffdiagonal, NULL); CHKERRV(ierr);
  
  // initialize name
  name = "D";
  ierr = PetscObjectSetName((PetscObject)matrixD, name.c_str()); CHKERRV(ierr);

  // --------
  // create matrix E
  Mat matrixE;
  ierr = MatCreate(mpiCommunicator, &matrixE); CHKERRV(ierr);

  // initialize size
  nRowsLocal = this->finiteElementMethodFat_.data().functionSpace()->nDofsLocalWithoutGhosts() - nSharedDofsLocal_;
  nColumnsLocal = this->finiteElementMethodDiffusionTotal_.data().functionSpace()->nDofsLocalWithoutGhosts();
  ierr = MatSetSizes(matrixE, nRowsLocal, nColumnsLocal, PETSC_DECIDE, PETSC_DECIDE); CHKERRV(ierr);
  ierr = MatSetType(matrixE, MATAIJ); CHKERRV(ierr);
  
  ierr = MatSeqAIJSetPreallocation(matrixE, nNonZerosDiagonal, NULL); CHKERRV(ierr);
  ierr = MatMPIAIJSetPreallocation(matrixE, nNonZerosDiagonal, NULL, nNonZerosOffdiagonal, NULL); CHKERRV(ierr);
  
  // initialize name
  name = "E";
  ierr = PetscObjectSetName((PetscObject)matrixE, name.c_str()); CHKERRV(ierr);

  // output sizes for debugging
  PetscInt nRowsMatrixB, nColumnsMatrixB;
  PetscInt nRowsMatrixC, nColumnsMatrixC;
  PetscInt nRowsMatrixD, nColumnsMatrixD;
  PetscInt nRowsMatrixE, nColumnsMatrixE;
  ierr = MatGetSize(matrixB, &nRowsMatrixB, &nColumnsMatrixB); CHKERRV(ierr);
  ierr = MatGetSize(matrixC, &nRowsMatrixC, &nColumnsMatrixC); CHKERRV(ierr);
  ierr = MatGetSize(matrixD, &nRowsMatrixD, &nColumnsMatrixD); CHKERRV(ierr);
  ierr = MatGetSize(matrixE, &nRowsMatrixE, &nColumnsMatrixE); CHKERRV(ierr);
  
  LOG(DEBUG) << "n shared nodes local: " << sharedNodes_.size() << ", n shared dofs local: " << nSharedDofsLocal_ << ", global: " << borderDofsGlobalFat_.size() 
    << ", create gamma matrices: matrixB: " << nRowsMatrixB << "x" << nColumnsMatrixB << ", matrixC: " << nRowsMatrixC << "x" << nColumnsMatrixC 
    << ", matrixD: " << nRowsMatrixD << "x" << nColumnsMatrixD << ", matrixE: " << nRowsMatrixE << "x" << nColumnsMatrixE;

  // compute entries
  setEntriesBorderMatrices(originalMatrixB, originalMatrixC, matrixB, matrixC, matrixD, matrixE);

  // store the matrices in the system matrix
  this->submatricesSystemMatrix_[(this->nCompartments_+0)*this->nColumnSubmatricesSystemMatrix_ + this->nCompartments_+0] = matrixB;
  this->submatricesSystemMatrix_[(this->nCompartments_+1)*this->nColumnSubmatricesSystemMatrix_ + this->nCompartments_+1] = matrixC;
  this->submatricesSystemMatrix_[(this->nCompartments_+0)*this->nColumnSubmatricesSystemMatrix_ + this->nCompartments_+1] = matrixD;
  this->submatricesSystemMatrix_[(this->nCompartments_+1)*this->nColumnSubmatricesSystemMatrix_ + this->nCompartments_+0] = matrixE;
}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusionMuscle,typename FiniteElementMethodDiffusionFat>
void MultidomainWithFatSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle,FiniteElementMethodDiffusionFat>::
setEntriesBorderMatrices(Mat originalMatrixB, Mat originalMatrixC, Mat matrixB, Mat matrixC, Mat matrixD, Mat matrixE)
{
  PetscErrorCode ierr;
  MPI_Comm mpiCommunicator = this->dataFat_.functionSpace()->meshPartition()->mpiCommunicator();
  PetscInt nColumnsGlobal = this->dataFat_.functionSpace()->nDofsGlobal();

#ifndef NDEBUG
  PetscInt nColumnsOriginalMatrixC = 0;
  ierr = MatGetSize(originalMatrixC, NULL, &nColumnsOriginalMatrixC); CHKERRV(ierr);
  assert(nColumnsOriginalMatrixC == nColumnsGlobal);

  PetscInt nRowsMatrixC = 0;
  PetscInt nColumnsMatrixC = 0;
  ierr = MatGetSize(matrixC, &nRowsMatrixC, &nColumnsMatrixC); CHKERRV(ierr);
  LOG(DEBUG) << "size of matrix C: " << nRowsMatrixC << "x" << nColumnsMatrixC;
#endif  

  // iterate over non-zero entries of the matrix originalMatrixC, this has to be done by creating the local sub matrix from which only the non-zero entries can be retrieved by MatGetRow
  PetscInt rowNoGlobalBegin = 0;
  PetscInt rowNoGlobalEnd = 0;
  ierr = MatGetOwnershipRange(originalMatrixC, &rowNoGlobalBegin, &rowNoGlobalEnd); CHKERRV(ierr);

  PetscInt nRowsLocal = rowNoGlobalEnd - rowNoGlobalBegin;

  // create index set indicating all rows of originalMatrixC that are stored locally, in global numbering
  IS indexSetRows[1];
  ISCreateStride(mpiCommunicator, nRowsLocal, rowNoGlobalBegin, 1, &indexSetRows[0]);

  // create index set indicating all columns of originalMatrixC that are stored locally (which are also all global columns), in global numbering
  IS indexSetColumns[1];
  ISCreateStride(mpiCommunicator, nColumnsGlobal, 0, 1, &indexSetColumns[0]);
  
  Mat *localSubMatrix;
  ierr = MatCreateSubMatrices(originalMatrixC, 1, indexSetRows, indexSetColumns, MAT_INITIAL_MATRIX, &localSubMatrix); CHKERRV(ierr);

  PetscInt nRowsLocalSubMatrix = 0;
  PetscInt nColumnsLocalSubMatrix = 0;
  ierr = MatGetSize(localSubMatrix[0], &nRowsLocalSubMatrix, &nColumnsLocalSubMatrix); CHKERRV(ierr);

  LOG(DEBUG) << "originalMatrixC has " << nRowsLocal << " local rows: [" << rowNoGlobalBegin << "," << rowNoGlobalEnd 
    << "], local submatrix is " << nRowsLocalSubMatrix << "x" << nColumnsLocalSubMatrix;

  // loop over rows of originalMatrixC
  for (PetscInt rowNoGlobal = rowNoGlobalBegin; rowNoGlobal < rowNoGlobalEnd; rowNoGlobal++)
  {
    PetscInt rowNoLocal = rowNoGlobal - rowNoGlobalBegin;
    
    // get non-zero values of current row
    PetscInt nNonzeroEntriesInRow;
    const PetscInt *columnIndices;
    const double *values;
    ierr = MatGetRow(localSubMatrix[0], rowNoLocal, &nNonzeroEntriesInRow, &columnIndices, &values); CHKERRV(ierr);

    bool currentRowDofIsBorder = borderDofsFat_.find(rowNoLocal) != borderDofsFat_.end();

    std::stringstream s;
    // loop over columns
    for (int columnIndex = 0; columnIndex < nNonzeroEntriesInRow; columnIndex++)
    {
      s << columnIndices[columnIndex] << " ";
    }
    VLOG(1) << "  row " << rowNoGlobal << " has " << nNonzeroEntriesInRow << " entries at " << s.str();

    // loop over columns
    for (int columnIndex = 0; columnIndex < nNonzeroEntriesInRow; columnIndex++)
    {
      PetscInt columnNoGlobal = columnIndices[columnIndex];
      double value = values[columnIndex];

      bool currentColumnDofIsBorder = borderDofsGlobalFat_.find(columnNoGlobal) != borderDofsGlobalFat_.end();

      VLOG(2) << "C[" << rowNoGlobal << "," << columnNoGlobal << "] = " << value;

      // if the current row is a shared dof
      if (currentRowDofIsBorder)
      {
        // if the current column is a shared dof
        if (currentColumnDofIsBorder)
        {
          PetscInt rowNoGlobalMuscle = getDofNoGlobalMuscleFromDofNoGlobalFat(rowNoGlobal);
          PetscInt columnNoGlobalMuscle = getDofNoGlobalMuscleFromDofNoGlobalFat(columnNoGlobal);

          // add the matrix entry of C to matrix B
          if (enableFatComputation_)
            ierr = MatSetValue(matrixB, rowNoGlobalMuscle, columnNoGlobalMuscle, value, ADD_VALUES); CHKERRV(ierr);
        }
        else 
        {
          // the current entry if from a shared row but not shared column
          PetscInt columnNoGlobalWithoutSharedDofs = getDofNoGlobalFatWithoutSharedDofs(columnNoGlobal);
          PetscInt rowNoGlobalMuscle = getDofNoGlobalMuscleFromDofNoGlobalFat(rowNoGlobal);

          // set the entry in matrix D
          ierr = MatSetValue(matrixD, rowNoGlobalMuscle, columnNoGlobalWithoutSharedDofs, value, INSERT_VALUES); CHKERRV(ierr);
        }
      }
      else 
      {
        // if the current row is not a shared dof
        
        // if the current column is a shared dof
        if (currentColumnDofIsBorder)
        {
          // the current entry if from a shared row but not shared column
          PetscInt rowNoGlobalWithoutSharedDofs = getDofNoGlobalFatWithoutSharedDofs(rowNoGlobal);
          PetscInt columnNoGlobalMuscle = getDofNoGlobalMuscleFromDofNoGlobalFat(columnNoGlobal);

          // set the entry in matrix E
          ierr = MatSetValue(matrixE, rowNoGlobalWithoutSharedDofs, columnNoGlobalMuscle, value, INSERT_VALUES); CHKERRV(ierr);
        }
        else 
        {
          // neither the row dof nor the column dof is shared, set the value in the new C matrix
          PetscInt rowNoGlobalWithoutSharedDofs = getDofNoGlobalFatWithoutSharedDofs(rowNoGlobal);
          PetscInt columnNoGlobalWithoutSharedDofs = getDofNoGlobalFatWithoutSharedDofs(columnNoGlobal);

          //LOG(DEBUG) << "fat full: (" << rowNoGlobal << "," << columnNoGlobal << ")/([" << rowNoGlobalBegin << "," << rowNoGlobalEnd << "]," << nColumnsGlobal 
          //  << ") without shared: (" << rowNoGlobalWithoutSharedDofs << "," << columnNoGlobalWithoutSharedDofs << ")/(" << nRowsMatrixC << "," << nColumnsMatrixC << ")";

          // set the entry in matrix C
          ierr = MatSetValue(matrixC, rowNoGlobalWithoutSharedDofs, columnNoGlobalWithoutSharedDofs, value, INSERT_VALUES); CHKERRV(ierr);
        }
      }
    }
  }

  // assembly of parallel matrices
  ierr = MatAssemblyBegin(matrixB, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyBegin(matrixC, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyBegin(matrixD, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyBegin(matrixE, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(matrixB, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(matrixC, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(matrixD, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(matrixE, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);

  LOG(DEBUG) << "setEntriesBorderMatrices";
  LOG(DEBUG) << "B: " << matrixB << ", C: " << matrixC << ", D: " << matrixD << ", E: " << matrixE;

  if (VLOG_IS_ON(1))
  {
    VLOG(1) << "B: " << matrixB << ", " << PetscUtility::getStringMatrix(matrixB);
    VLOG(1) << "C: " << matrixC << ", " << PetscUtility::getStringMatrix(matrixC);
    VLOG(1) << "D: " << matrixD << ", " << PetscUtility::getStringMatrix(matrixD);
    VLOG(1) << "E: " << matrixE << ", " << PetscUtility::getStringMatrix(matrixE);
  }
}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusionMuscle,typename FiniteElementMethodDiffusionFat>
void MultidomainWithFatSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle,FiniteElementMethodDiffusionFat>::
updateBorderMatrices()
{
  // this method computes the entries again, after the geometry has updated by the contraction

  Mat originalMatrixB = this->finiteElementMethodDiffusionTotal_.data().stiffnessMatrix()->valuesGlobal();
  Mat originalMatrixC = this->finiteElementMethodFat_.data().stiffnessMatrix()->valuesGlobal();

  // store the matrices in the system matrix
  Mat &matrixB = this->submatricesSystemMatrix_[(this->nCompartments_+0)*this->nColumnSubmatricesSystemMatrix_ + this->nCompartments_+0];
  Mat &matrixC = this->submatricesSystemMatrix_[(this->nCompartments_+1)*this->nColumnSubmatricesSystemMatrix_ + this->nCompartments_+1];
  Mat &matrixD = this->submatricesSystemMatrix_[(this->nCompartments_+0)*this->nColumnSubmatricesSystemMatrix_ + this->nCompartments_+1];
  Mat &matrixE = this->submatricesSystemMatrix_[(this->nCompartments_+1)*this->nColumnSubmatricesSystemMatrix_ + this->nCompartments_+0];

  // compute entries
  setEntriesBorderMatrices(originalMatrixB, originalMatrixC, matrixB, matrixC, matrixD, matrixE);
}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusionMuscle,typename FiniteElementMethodDiffusionFat>
void MultidomainWithFatSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle,FiniteElementMethodDiffusionFat>::
copyPhiBToSolution()
{
  // get entries from: dataFat_.extraCellularPotentialFat()->valuesGlobal(0)
  // set entries in:  this->subvectorsSolution_[this->nCompartments_+1]
  // discard the entries for shared dofs

  Vec phiB = dataFat_.extraCellularPotentialFat()->valuesGlobal(0);
  Vec solution = this->subvectorsSolution_[this->nCompartments_+1];

  // get the local range that is stored on this rank for phiB
  PetscInt dofNoGlobalBegin, dofNoGlobalEnd;
  PetscErrorCode ierr;
  ierr = VecGetOwnershipRange(phiB, &dofNoGlobalBegin, &dofNoGlobalEnd); CHKERRV(ierr);

  PetscInt nEntriesLocal = dofNoGlobalEnd - dofNoGlobalBegin;

  // get all values in phiB
  std::vector<PetscInt> indices(nEntriesLocal);
  std::iota(indices.begin(), indices.end(), dofNoGlobalBegin);
  std::vector<double> values(nEntriesLocal);
  ierr = VecGetValues(phiB, nEntriesLocal, indices.data(), values.data()); CHKERRV(ierr);

  // loop over local entries of fat
  for (PetscInt dofNoGlobal = dofNoGlobalBegin; dofNoGlobal < dofNoGlobalEnd; dofNoGlobal++)
  {
    PetscInt dofNoLocalFat = dofNoGlobal - dofNoGlobalBegin;

    // if dof is not shared
    if (borderDofsFat_.find(dofNoLocalFat) == borderDofsFat_.end())
    {
      // set the value in the solution
      double value = values[dofNoLocalFat];
      PetscInt dofNoGlobalWithoutSharedDofs = getDofNoGlobalFatWithoutSharedDofs(dofNoGlobal);
      ierr = VecSetValue(solution, dofNoGlobalWithoutSharedDofs, value, INSERT_VALUES); CHKERRV(ierr);
    }
  }

  // global assembly of solution Vec
  ierr = VecAssemblyBegin(solution); CHKERRV(ierr);
  ierr = VecAssemblyEnd(solution); CHKERRV(ierr);

  LOG(DEBUG) << "copyPhiBToSolution: get from phi_b and set in solution";
  LOG(DEBUG) << "borderDofsFat_: " << borderDofsFat_;
  LOG(DEBUG) << "phi_b: " << PetscUtility::getStringVector(phiB);
  LOG(DEBUG) << "solution: " << PetscUtility::getStringVector(solution);
}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusionMuscle,typename FiniteElementMethodDiffusionFat>
void MultidomainWithFatSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle,FiniteElementMethodDiffusionFat>::
copySolutionToPhiB()
{
  // get entries from: this->subvectorsSolution_[this->nCompartments_+1]  (phi_b) (does not contain values for shared dofs)
  // and:              this->subvectorsSolution_[this->nCompartments_]    (phi_e) (only get the shared dof values from here)
  // set entries in:  dataFat_.extraCellularPotentialFat()->valuesGlobal(0)

  Vec phiB = this->subvectorsSolution_[this->nCompartments_+1];
  Vec phiE = this->subvectorsSolution_[this->nCompartments_];
  Vec result = dataFat_.extraCellularPotentialFat()->valuesGlobal(0);

  // get the local range that of phiB is stored on this rank for result
  PetscInt dofNoGlobalBeginPhiB, dofNoGlobalEndPhiB;
  PetscErrorCode ierr;
  ierr = VecGetOwnershipRange(phiB, &dofNoGlobalBeginPhiB, &dofNoGlobalEndPhiB); CHKERRV(ierr);

  PetscInt nEntriesLocal = dofNoGlobalEndPhiB - dofNoGlobalBeginPhiB;

  // get all values in phiB
  std::vector<PetscInt> indicesPhiB(nEntriesLocal);
  std::iota(indicesPhiB.begin(), indicesPhiB.end(), dofNoGlobalBeginPhiB);
  std::vector<double> valuesPhiB(nEntriesLocal);
  
  ierr = VecGetValues(phiB, nEntriesLocal, indicesPhiB.data(), valuesPhiB.data()); CHKERRV(ierr);

  // get the local range of phiE that is stored on this rank for result
  PetscInt dofNoGlobalBeginPhiE, dofNoGlobalEndPhiE;
  ierr = VecGetOwnershipRange(phiE, &dofNoGlobalBeginPhiE, &dofNoGlobalEndPhiE); CHKERRV(ierr);
  nEntriesLocal = dofNoGlobalEndPhiE - dofNoGlobalBeginPhiE;

  // get all values in phiE
  std::vector<PetscInt> indicesPhiE(nEntriesLocal);
  std::iota(indicesPhiE.begin(), indicesPhiE.end(), dofNoGlobalBeginPhiE);
  std::vector<double> valuesPhiE(nEntriesLocal);

  ierr = VecGetValues(phiE, nEntriesLocal, indicesPhiE.data(), valuesPhiE.data()); CHKERRV(ierr);

  PetscInt dofNoGlobalBeginFat, dofNoGlobalEndFat;
  ierr = VecGetOwnershipRange(result, &dofNoGlobalBeginFat, &dofNoGlobalEndFat); CHKERRV(ierr);

  LOG(DEBUG) << "phiB: [" << dofNoGlobalBeginPhiB << "," << dofNoGlobalEndPhiB << "]";
  LOG(DEBUG) << "phiE: [" << dofNoGlobalBeginPhiE << "," << dofNoGlobalEndPhiE << "]";
  LOG(DEBUG) << "result: [" << dofNoGlobalBeginFat << "," << dofNoGlobalEndFat << "]";

  // loop over local entries of result which is on the fat mesh
  for (PetscInt dofNoGlobalFat = dofNoGlobalBeginFat; dofNoGlobalFat < dofNoGlobalEndFat; dofNoGlobalFat++)
  {
    PetscInt dofNoLocalFat = dofNoGlobalFat - dofNoGlobalBeginFat;

    // if dof is shared, use value from phiE
    if (borderDofsFat_.find(dofNoLocalFat) != borderDofsFat_.end())
    {
      // set the value in the result
      PetscInt dofNoGlobalMuscle = getDofNoGlobalMuscleFromDofNoGlobalFat(dofNoGlobalFat);
      PetscInt dofNoLocalMuscle = dofNoGlobalMuscle - dofNoGlobalBeginPhiE;
      double value = valuesPhiE[dofNoLocalMuscle];
      ierr = VecSetValue(result, dofNoGlobalFat, value, INSERT_VALUES); CHKERRV(ierr);
    }
    else 
    {
      // if dof is not shared, use value from phiB
      PetscInt dofNoGlobalFatWithoutSharedDofs = getDofNoGlobalFatWithoutSharedDofs(dofNoGlobalFat);
      PetscInt dofNoLocalFat = dofNoGlobalFatWithoutSharedDofs - dofNoGlobalBeginPhiB;
      double value = valuesPhiB[dofNoLocalFat];
      ierr = VecSetValue(result, dofNoGlobalFat, value, INSERT_VALUES); CHKERRV(ierr);
    }
  }

  // global assembly of result Vec
  ierr = VecAssemblyBegin(result); CHKERRV(ierr);
  ierr = VecAssemblyEnd(result); CHKERRV(ierr);

  LOG(DEBUG) << "copySolutionToPhiB: from phi_b and phi_e get result = whole phi_b";
  LOG(DEBUG) << "borderDofsFat_: " << borderDofsFat_;
  LOG(DEBUG) << "phi_b (without shared dofs): " << PetscUtility::getStringVector(phiB);
  LOG(DEBUG) << "phi_e (with shared dofs): " << PetscUtility::getStringVector(phiE);
  LOG(DEBUG) << "result: " << PetscUtility::getStringVector(result);
}

} // namespace TimeSteppingScheme
