#include "specialized_solver/multidomain_solver/multidomain_with_fat_solver.h"

#include <Python.h>  // has to be the first included header

namespace TimeSteppingScheme
{

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusionMuscle,typename FiniteElementMethodDiffusionFat>
MultidomainWithFatSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle,FiniteElementMethodDiffusionFat>::
MultidomainWithFatSolver(DihuContext context) :
  MultidomainSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle>(context),
  dataFat_(this->context_),
  finiteElementMethodFat_(this->context_["Fat"])
{
}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusionMuscle,typename FiniteElementMethodDiffusionFat>
void MultidomainWithFatSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle,FiniteElementMethodDiffusionFat>::
initialize()
{
  MultidomainSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle>::initializeObjects();

  // indicate in solverStructureVisualizer that now a child solver will be initialized
  DihuContext::solverStructureVisualizer()->beginChild("Fat");

  LOG(DEBUG) << "initialize fat FEM";
  // initialize the potential flow finite element method, this also creates the function space
  finiteElementMethodFat_.initialize();

  // indicate in solverStructureVisualizer that the child solver initialization is done
  DihuContext::solverStructureVisualizer()->endChild();

  LOG(DEBUG) << "function space \"" << finiteElementMethodFat_.functionSpace()->meshName() << "\".";

  // initialize the data object
  dataFat_.setFunctionSpace(finiteElementMethodFat_.functionSpace());
  dataFat_.setDataMultidomain(
    std::make_shared<typename MultidomainSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle>::Data>(
    this->dataMultidomain_));
  
  dataFat_.initialize();

  // get θ value for Crank-Nicolson scheme
  theta_ = this->specificSettings_.getOptionDouble("theta", 0.5, PythonUtility::Positive);
  useLumpedMassMatrix_ = this->specificSettings_.getOptionBool("useLumpedMassMatrix", true);
  enableFatComputation_ = this->specificSettings_.getOptionBool("enableFatComputation", true);
  if (!enableFatComputation_)
  {
    LOG(WARNING) << this->specificSettings_ << "[\"enableFatComputation\"] is set to false. This will disable the fat layer computation.";
  }

  // initialize sharedNodes_ i.e. the border nodes that are shared between muscle and fat mesh
  findSharedNodesBetweenMuscleAndFat();

  // system to be solved (here for nCompartments_=3):
  //
  // [A^1_Vm,Vm   |            |             | B^1_Vm,phie |             ]   [ V^1_m^(i+1) ]    [b^1^(i)]  <── n rows: number of dofs in muscle mesh
  // [            | A^2_Vm,Vm  |             | B^2_Vm,phie |             ]   [ V^2_m^(i+1) ]    [b^2^(i)]  <── n rows: number of dofs in muscle mesh
  // [   ...      |            | A^M_Vm,Vm   | B^M_Vm,phie |             ] * [ V^M_m^(i+1) ] =  [b^M^(i)]  <── n rows: number of dofs in muscle mesh
  // [B^1_phie,Vm |B^2_phie,Vm | B^M_phie,Vm | B_phie,phie |      D      ]   [ phi_e^(i+1) ]    [ 0     ]  <── n rows: number of dofs in muscle mesh
  // [            |            |             |     E       | C_phib,phib ]   [ phi_b^(i+1) ]    [ 0     ]  <── n rows: number of dofs in fat mesh minus number of shared dofs
  //     ^              ^          ^               ^              ^
  //     |              |          |               |              |
  //  -----n------cols: n dofs in muscle mesh--------------     n cols: n dofs in fat mesh minus number of shared dofs
  //
  // see multidomain.pdf for the definition of the submatrices, B_phie,phie is also called B and C_phib,phib is called C

  // initialize rhs and solution vector, the entries up to this->nCompartments_ were already set by the MultidomainSolver
  this->subvectorsRightHandSide_.resize(this->nCompartments_+2);
  this->subvectorsSolution_.resize(this->nCompartments_+2);

  // initialize all submatrices for system matrix, except B,C,D,E
  this->timeStepWidthOfSystemMatrix_ = this->timeStepWidth_;
  setSystemMatrixSubmatrices(this->timeStepWidthOfSystemMatrix_);

  // initialize matices B,C,D,E in submatrices for system matrix,
  // also set the last zero entry of the rhs and the entry for phi_b^(i+1) in the solution vector
  initializeBorderVariables();

  // from the initialize submatrices create the actual system matrix, this->nestedSystemMatrix_ 
  // as nested Petsc Mat and also the single Mat, this->singleSystemMatrix_ 
  this->createSystemMatrixFromSubmatrices();

  // set the nullspace of the matrix
  // as we have Neumann boundary conditions, constant functions are in the nullspace of the matrix
  MatNullSpace nullSpace;
  PetscErrorCode ierr;
  ierr = MatNullSpaceCreate(data().functionSpace()->meshPartition()->mpiCommunicator(), PETSC_TRUE, 0, PETSC_NULL, &nullSpace); CHKERRV(ierr);
  ierr = MatSetNullSpace(this->singleSystemMatrix_, nullSpace); CHKERRV(ierr);
  ierr = MatSetNearNullSpace(this->singleSystemMatrix_, nullSpace); CHKERRV(ierr); // for multigrid methods
  //ierr = MatNullSpaceDestroy(&nullSpace); CHKERRV(ierr);

  // set matrix used for linear solver and preconditioner to ksp context
  assert(this->linearSolver_->ksp());
  ierr = KSPSetOperators(*this->linearSolver_->ksp(), this->singleSystemMatrix_, this->singlePreconditionerMatrix_); CHKERRV(ierr);

  if (this->alternativeLinearSolver_)
    ierr = KSPSetOperators(*this->alternativeLinearSolver_->ksp(), this->singleSystemMatrix_, this->singlePreconditionerMatrix_); CHKERRV(ierr);

  // set block information in preconditioner for block jacobi and node positions for MG preconditioners
  setInformationToPreconditioner();

  // create temporary vector which is needed in computation of rhs, b1_ was created by setSystemMatrixSubmatrices()
  ierr = MatCreateVecs(b1_[0], &temporary_, NULL); CHKERRV(ierr);

  // initialize subvectors for rhs
  // entry nCompartments was already set to zero by MultidomainSolver, entry nCompartments+1 was already set by initializeBorderVariables()
  for (int k = 0; k < this->nCompartments_; k++)
  {    
    // initialize top part of rhs
    ierr = MatCreateVecs(b1_[k], &this->subvectorsRightHandSide_[k], NULL); CHKERRV(ierr);
  }

  // set values for phi_e
  this->subvectorsRightHandSide_[this->nCompartments_] = this->dataMultidomain_.zero()->valuesGlobal();
  // subvectorsRightHandSide_[nCompartments_+1] has been set by initializeBorderVariables()

  // set vectors of Vm in compartments in subvectorsSolution_
  for (int k = 0; k < this->nCompartments_; k++)
  {
    this->subvectorsSolution_[k] = this->dataMultidomain_.transmembranePotentialSolution(k)->valuesGlobal(0); // this is for V_mk^(i+1)
    ierr = VecDuplicate(this->subvectorsSolution_[k], &this->subvectorsRightHandSide_[k]); CHKERRV(ierr);
    ierr = VecZeroEntries(this->subvectorsRightHandSide_[k]); CHKERRV(ierr);
  }
  
  // set values for phi_e
  this->subvectorsSolution_[this->nCompartments_] = this->dataMultidomain_.extraCellularPotential()->valuesGlobal();
  ierr = VecZeroEntries(this->subvectorsSolution_[this->nCompartments_]); CHKERRV(ierr);
  // subvectorsSolution_[nCompartments_+1] has been set by initializeBorderVariables()

  // create the nested Petsc Vec's
  LOG(DEBUG) << "create nested vector";
  ierr = VecCreateNest(this->rankSubset_->mpiCommunicator(), this->nCompartments_+2, NULL,
                       this->subvectorsRightHandSide_.data(), &this->nestedRightHandSide_); CHKERRV(ierr);

  ierr = VecCreateNest(this->rankSubset_->mpiCommunicator(), this->nCompartments_+2, NULL, 
                       this->subvectorsSolution_.data(), &this->nestedSolution_); CHKERRV(ierr);

  // write initial meshes
  callOutputWriter(0, 0.0, 0);
}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusionMuscle,typename FiniteElementMethodDiffusionFat>
void MultidomainWithFatSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle,FiniteElementMethodDiffusionFat>::
setSystemMatrixSubmatrices(double timeStepWidth)
{
  LOG(INFO) << "dt of system matrix: " << timeStepWidth;
  assert(this->finiteElementMethodDiffusionTotal_.data().stiffnessMatrix());
  assert(this->finiteElementMethodDiffusion_.data().stiffnessMatrix());
  assert(this->finiteElementMethodFat_.data().stiffnessMatrix());
  PetscErrorCode ierr;

  // initialize number of submatrix rows in the system matrix
  this->nColumnSubmatricesSystemMatrix_ = this->nCompartments_+1+1;

  this->submatricesSystemMatrix_.resize(MathUtility::sqr(this->nColumnSubmatricesSystemMatrix_),NULL);

  // system to be solved:
  //
  // [A^1_Vm,Vm   |            |             | B^1_Vm,phie |             ]   [ V^1_m^(i+1) ]    [b^1^(i)]
  // [            | A^2_Vm,Vm  |             | B^2_Vm,phie |             ]   [ V^2_m^(i+1) ]    [b^2^(i)]
  // [   ...      |            | A^M_Vm,Vm   | B^M_Vm,phie |             ] * [ V^M_m^(i+1) ] =  [b^M^(i)]
  // [B^1_phie,Vm |B^2_phie,Vm | B^M_phie,Vm | B_phie,phie |      D      ]   [ phi_e^(i+1) ]    [ 0     ]
  // [            |            |             |     E       | C_phib,phib ]   [ phi_b^(i+1) ]    [ 0     ]

  LOG(TRACE) << "setSystemMatrix";

  // fill this->submatricesSystemMatrix_, empty submatrices may be NULL
  // stiffnessMatrix and inverse lumped mass matrix without prefactor
  Mat stiffnessMatrix = this->finiteElementMethodDiffusion_.data().stiffnessMatrix()->valuesGlobal();
  Mat massMatrix = this->finiteElementMethodDiffusion_.data().massMatrix()->valuesGlobal();

  Mat minusDtMInv;   // matrix -dt*M^-1 with which to scale the matrix equations of the compartments if option useLumpedMassMatrix_ is set
  if (useLumpedMassMatrix_)
  {
    Mat inverseLumpedMassMatrix = this->finiteElementMethodDiffusion_.data().inverseLumpedMassMatrix()->valuesGlobal();
    
    // compute minusDtMInv = -dt*M^-1
    ierr = MatConvert(inverseLumpedMassMatrix, MATSAME, MAT_INITIAL_MATRIX, &minusDtMInv); CHKERRV(ierr);
    ierr = MatScale(minusDtMInv, -timeStepWidth); CHKERRV(ierr);
  }

  // set all submatrices
  for (int k = 0; k < this->nCompartments_; k++)
  {
    // right column matrix
    double prefactor = theta_ / (this->am_[k]*this->cm_[k]);

    VLOG(2) << "k=" << k << ", am: " << this->am_[k] << ", cm: " << this->cm_[k] << ", prefactor: " << prefactor;

    // matrix B on right column
    // create matrix as theta/(Am*Cm)*K
    Mat matrixOnRightColumn;
    ierr = MatConvert(stiffnessMatrix, MATSAME, MAT_INITIAL_MATRIX, &matrixOnRightColumn); CHKERRV(ierr);

    // scale matrix on right column with prefactor
    ierr = MatScale(matrixOnRightColumn, prefactor); CHKERRV(ierr);

    if (useLumpedMassMatrix_)
    {
      // in this formulation the matrix B is B = -dt*theta/(Am*Cm)*M^-1*K
      Mat result;
      ierr = MatMatMult(minusDtMInv, matrixOnRightColumn, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result); CHKERRV(ierr);
      matrixOnRightColumn = result;
    }

    // set on right column of the system matrix
    this->submatricesSystemMatrix_[k*this->nColumnSubmatricesSystemMatrix_ + this->nCompartments_] = matrixOnRightColumn;

    // ---
    // matrix on diagonal, A
    // copy right block matrix also to diagonal matrix
    Mat matrixOnDiagonalBlock;
    ierr = MatConvert(matrixOnRightColumn, MATSAME, MAT_INITIAL_MATRIX, &matrixOnDiagonalBlock); CHKERRV(ierr);

    if (useLumpedMassMatrix_)
    {
      // in this formulation the matrix A is A = -dt*theta/(Amk*Cmk)*M^-1*K + I, with B = -dt*theta/(Amk*Cmk)*M^-1*K this becomes A = B + I
      // add identity
      ierr = MatShift(matrixOnDiagonalBlock, 1); CHKERRV(ierr);
    }
    else 
    {
      // in this formulation the matrix A is A = theta/(Amk*Cmk)*K - 1/dt*M, with B = theta/(Am*Cm)*K this becomes A = B - 1/dt*M

      // add scaled mass matrix, -1/dt*M,  AXPY: Y = a*X + Y  MatAXPY(Y,a,X,SAME_NONZERO_PATTERN)
      prefactor = -1/timeStepWidth;
      ierr = MatAXPY(matrixOnDiagonalBlock, prefactor, massMatrix, SAME_NONZERO_PATTERN); CHKERRV(ierr);
    }

    // set on diagonal
    this->submatricesSystemMatrix_[k*this->nColumnSubmatricesSystemMatrix_ + k] = matrixOnDiagonalBlock;

    // ---
    // bottom row matrices, B
    // stiffnessMatrixWithPrefactor is f_k*K
    Mat stiffnessMatrixWithPrefactor = this->finiteElementMethodDiffusionCompartment_[k].data().stiffnessMatrix()->valuesGlobal();

    // create matrix as copy of stiffnessMatrix
    Mat matrixOnBottomRow;
    ierr = MatConvert(stiffnessMatrixWithPrefactor, MATSAME, MAT_INITIAL_MATRIX, &matrixOnBottomRow); CHKERRV(ierr);

    // set on bottom row of the system matrix
    this->submatricesSystemMatrix_[this->nCompartments_*this->nColumnSubmatricesSystemMatrix_ + k] = matrixOnBottomRow;
  }

  // ---
  // matrices for rhs
  b1_.resize(this->nCompartments_);
  b2_.resize(this->nCompartments_);

  // initialize the matrices b1, b2 to compute the rhs
  // set all this->submatricesSystemMatrix_
  for (int k = 0; k < this->nCompartments_; k++)
  {
    // set b1_ = (θ-1)*1/(Am^k*Cm^k)*K_sigmai^k
    ierr = MatConvert(stiffnessMatrix, MATSAME, MAT_INITIAL_MATRIX, &b1_[k]); CHKERRV(ierr);
    
    double prefactor = (theta_ - 1) / (this->am_[k]*this->cm_[k]);
    ierr = MatScale(b1_[k], prefactor); CHKERRV(ierr);

    if (useLumpedMassMatrix_)
    {
      // in this formulation we have b1_[k] = -dt*(θ-1)/(Am^k*Cm^k)*M^{-1}*K_sigmai^k + I
      Mat result;
      ierr = MatMatMult(minusDtMInv, b1_[k], MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result); CHKERRV(ierr);
      b1_[k] = result;
      
      // add identity
      ierr = MatShift(b1_[k], 1); CHKERRV(ierr);
    }
    else 
    {
      // add scaled mass matrix, -1/dt*M,  AXPY: Y = a*X + Y  MatAXPY(Y,a,X,SAME_NONZERO_PATTERN)
      prefactor = -1/timeStepWidth;
      ierr = MatAXPY(b1_[k], prefactor, massMatrix, SAME_NONZERO_PATTERN); CHKERRV(ierr);
    }

    // set b2_ = (θ-1)*K_sigmai^k
    ierr = MatConvert(stiffnessMatrix, MATSAME, MAT_INITIAL_MATRIX, &b2_[k]); CHKERRV(ierr);
    
    if (useLumpedMassMatrix_)
    {
      // in this formulation we have b2_[k] = -dt*(θ-1)*M^{-1}*K_sigmai^k
      Mat result;
      ierr = MatMatMult(minusDtMInv, b2_[k], MAT_INITIAL_MATRIX, PETSC_DEFAULT, &result); CHKERRV(ierr);
      b2_[k] = result;
    }

    prefactor = theta_ - 1;
    ierr = MatScale(b2_[k], prefactor); CHKERRV(ierr);
  }
}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusionMuscle,typename FiniteElementMethodDiffusionFat>
void MultidomainWithFatSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle,FiniteElementMethodDiffusionFat>::
setInformationToPreconditioner()
{
  // set block information in preconditioner, if it is of type block jacobi
  PetscErrorCode ierr;
  PC pc;
  ierr = KSPGetPC(*this->linearSolver_->ksp(), &pc); CHKERRV(ierr);
  
  // set block information for block jacobi preconditioner
  // check, if block jacobi preconditioner is selected
  PetscBool useBlockJacobiPreconditioner;
  PetscObjectTypeCompare((PetscObject)pc, PCBJACOBI, &useBlockJacobiPreconditioner);

  if (useBlockJacobiPreconditioner)
  {
    // smaller blocks
#if 1
    int nRanks = this->dataMultidomain_.functionSpace()->meshPartition()->rankSubset()->size();
    PetscInt nMatrixBlocks = this->nColumnSubmatricesSystemMatrix_*nRanks;
    
    std::shared_ptr<Partition::MeshPartition<FunctionSpace>> meshPartitionMuscle = this->dataMultidomain_.functionSpace()->meshPartition();
    std::shared_ptr<Partition::MeshPartition<typename FiniteElementMethodDiffusionFat::FunctionSpace>> meshPartitionFat = this->finiteElementMethodFat_.data().functionSpace()->meshPartition();

    // set sizes of all blocks to the number of dofs in the muscle domain
    std::vector<PetscInt> lengthsOfBlocks(nMatrixBlocks, this->dataMultidomain_.functionSpace()->nDofsLocalWithoutGhosts());
    
    // loop over ranks
    for (int rankNo = 0; rankNo < nRanks; rankNo++)
    {
      // loop over matrix blocks of muscle mesh, not the fat mesh
      for (int blockIndex = 0; blockIndex < this->nColumnSubmatricesSystemMatrix_-1; blockIndex++)
      {
        // get local size of block on the current rank
        // for muscle mesh block
        // determine number of local nodes (= dofs) for rank rankNo
        int nNodesLocalWithoutGhosts = 1;
        for (int coordinateDirection = 0; coordinateDirection < 3; coordinateDirection++)
        {
          int partitionIndex = meshPartitionMuscle->convertRankNoToPartitionIndex(coordinateDirection, rankNo);
          nNodesLocalWithoutGhosts *= meshPartitionMuscle->nNodesLocalWithoutGhosts(coordinateDirection, partitionIndex);
          VLOG(1) << "  block " << blockIndex << " rank " << rankNo << " dim " << coordinateDirection << ": *" << meshPartitionMuscle->nNodesLocalWithoutGhosts(coordinateDirection, rankNo) << " -> " << nNodesLocalWithoutGhosts;
        }

        VLOG(1) << "-> lengthsOfBlocks[" << rankNo*this->nColumnSubmatricesSystemMatrix_ + blockIndex << "] = " << nNodesLocalWithoutGhosts;
        lengthsOfBlocks[rankNo*this->nColumnSubmatricesSystemMatrix_ + blockIndex] = nNodesLocalWithoutGhosts;
      }
    }

    // set entries for block of fat meh
    Mat matrixC = this->submatricesSystemMatrix_[MathUtility::sqr(this->nColumnSubmatricesSystemMatrix_)-1];
    const PetscInt *matrixCOwnershipRanges;
    ierr = MatGetOwnershipRanges(matrixC, &matrixCOwnershipRanges); CHKERRV(ierr);

    // loop over ranks
    for (int rankNo = 0; rankNo < nRanks; rankNo++)
    {
      int nRowsLocalCurrentRank = matrixCOwnershipRanges[rankNo+1] - matrixCOwnershipRanges[rankNo];
      lengthsOfBlocks[rankNo*this->nColumnSubmatricesSystemMatrix_ + this->nColumnSubmatricesSystemMatrix_-1] = nRowsLocalCurrentRank;
    }
    

    // assert that size matches global matrix size
    PetscInt nRowsGlobal, nColumnsGlobal;
    ierr = MatGetSize(this->singleSystemMatrix_, &nRowsGlobal, &nColumnsGlobal); CHKERRV(ierr);
    PetscInt size = 0;
    for (PetscInt blockSize : lengthsOfBlocks)
      size += blockSize;

    LOG(DEBUG) << "block jacobi preconditioner, lengthsOfBlocks: " << lengthsOfBlocks << ", system matrix size: " << nRowsGlobal << "x" << nColumnsGlobal;

    if (size != nRowsGlobal || size != nColumnsGlobal)
    {
      LOG(FATAL) << "Block lengths of block jacobi preconditioner do not sum up to the system matrix size. Sum of block lengths: " << size << ", dimension of system matrix: "
        << nRowsGlobal << "x" << nColumnsGlobal;
    }
    assert(size == nRowsGlobal);
    assert(size == nColumnsGlobal);

    // PCBJacobiSetTotalBlocks(PC pc, PetscInt nBlocks, const PetscInt lengthsOfBlocks[])
    ierr = PCBJacobiSetTotalBlocks(pc, this->nColumnSubmatricesSystemMatrix_, lengthsOfBlocks.data()); CHKERRV(ierr);

#else
    // big blocks
    PetscInt nBlocks = this->nColumnSubmatricesSystemMatrix_;

    // set sizes of all blocks to the number of dofs in the muscle domain
    std::vector<PetscInt> lengthsOfBlocks(nBlocks, this->dataMultidomain_.functionSpace()->nDofsGlobal());
    
    // set last block size
    PetscInt nRowsGlobalLastSubMatrix;
    ierr = MatGetSize(this->submatricesSystemMatrix_[MathUtility::sqr(this->nColumnSubmatricesSystemMatrix_)-1], &nRowsGlobalLastSubMatrix, NULL); CHKERRV(ierr);
    lengthsOfBlocks[nBlocks-1] = nRowsGlobalLastSubMatrix;

    // assert that size matches global matrix size
    PetscInt nRowsGlobal, nColumnsGlobal;
    ierr = MatGetSize(this->singleSystemMatrix_, &nRowsGlobal, &nColumnsGlobal); CHKERRV(ierr);
    PetscInt size = 0;
    for (PetscInt blockSize : lengthsOfBlocks)
      size += blockSize;

    LOG(INFO) << "block jacobi preconditioner, lengthsOfBlocks: " << lengthsOfBlocks << ", system matrix size: " << nRowsGlobal << "x" << nColumnsGlobal;

    if (size != nRowsGlobal || size != nColumnsGlobal)
    {
      LOG(FATAL) << "Block lengths of block jacobi preconditioner do not sum up to the system matrix size.";
    }
    assert(size == nRowsGlobal);
    assert(size == nColumnsGlobal);

    // PCBJacobiSetTotalBlocks(PC pc, PetscInt nBlocks, const PetscInt lengthsOfBlocks[])
    ierr = PCBJacobiSetTotalBlocks(pc, this->nColumnSubmatricesSystemMatrix_, lengthsOfBlocks.data()); CHKERRV(ierr);
#endif
  }
  
  // set node positions

  // set the local node positions for the preconditioner
  int nNodesLocalMuscle = this->dataMultidomain_.functionSpace()->nNodesLocalWithoutGhosts();
  int nNodesLocalFat = this->dataFat_.functionSpace()->nNodesLocalWithoutGhosts() - nSharedDofsLocal_;
  int nNodesLocalBlock = nNodesLocalMuscle + nNodesLocalFat;

  std::vector<double> nodePositionCoordinatesForPreconditioner;
  nodePositionCoordinatesForPreconditioner.reserve(3*nNodesLocalBlock);

  // loop over muscle nodes and add their node positions
  for (dof_no_t dofNoLocalMuscle = 0; dofNoLocalMuscle < this->dataMultidomain_.functionSpace()->nDofsLocalWithoutGhosts(); dofNoLocalMuscle++)
  {
    Vec3 nodePosition = this->dataMultidomain_.functionSpace()->getGeometry(dofNoLocalMuscle);
  
    // add the coordinates
    for (int i = 0; i < 3; i++)
      nodePositionCoordinatesForPreconditioner.push_back(nodePosition[i]);
  }
  
  // loop over fat nodes and add their node positions
  for (dof_no_t dofNoLocalFat = 0; dofNoLocalFat < this->dataFat_.functionSpace()->nDofsLocalWithoutGhosts(); dofNoLocalFat++)
  {
    // if current fat dof is not shared
    if (borderDofsFat_.find(dofNoLocalFat) == borderDofsFat_.end())
    {
      Vec3 nodePosition = this->dataFat_.functionSpace()->getGeometry(dofNoLocalFat);
      
      // add the coordinates
      for (int i = 0; i < 3; i++)
        nodePositionCoordinatesForPreconditioner.push_back(nodePosition[i]);
    }
  }
  assert(nodePositionCoordinatesForPreconditioner.size() == 3*nNodesLocalBlock);

  LOG(DEBUG) << "set coordinates to preconditioner, " << nNodesLocalBlock << " node coordinates";

  ierr = PCSetCoordinates(pc, 3, nNodesLocalBlock, nodePositionCoordinatesForPreconditioner.data()); CHKERRV(ierr);

  // initialize preconditioner of alternative linear solver
  if (this->alternativeLinearSolver_)
  {
    PC pc;
    ierr = KSPGetPC(*this->alternativeLinearSolver_->ksp(), &pc); CHKERRV(ierr);
    ierr = PCSetCoordinates(pc, 3, nNodesLocalBlock, nodePositionCoordinatesForPreconditioner.data()); CHKERRV(ierr);
  }

  if (useBlockJacobiPreconditioner)
  {
    // for block jacobi set sub solvers
    // parse solver type of sub solver
    std::string subSolverType = this->specificSettings_.getOptionString("subSolverType", "none");
    std::string subPreconditionerType = this->specificSettings_.getOptionString("subPreconditionerType", "none");
    
    KSPType subKspType;
    PCType subPcType;
    Solver::Linear::parseSolverTypes(subSolverType, subPreconditionerType, subKspType, subPcType);
      
    // create internal data structures of ksp and pc, if not already done
    ierr = PCSetUp(pc); CHKERRV(ierr);

    // get sub ksp object
    KSP *subKspPointer;
    PetscInt nBlocksLocal;
    PetscInt firstBlockNoGlobal;
    ierr = PCBJacobiGetSubKSP(pc, &nBlocksLocal, &firstBlockNoGlobal, &subKspPointer); CHKERRV(ierr);
        
    LOG(INFO) << "using block jacobi, nBlocksLocal: " << nBlocksLocal << ", firstBlockNoGlobal: " << firstBlockNoGlobal 
      << ", subSolverType: " << subSolverType << ", subPreconditionerType: " << subPreconditionerType;
        
    for (int subKspLocalNo = 0; subKspLocalNo < nBlocksLocal; subKspLocalNo++)
    {
      KSP subKsp = subKspPointer[subKspLocalNo];
      
      // set solver type
      ierr = KSPSetType(subKsp, subKspType); CHKERRV(ierr);

      // set options from command line, this overrides the python config
      ierr = KSPSetFromOptions(subKsp); CHKERRV(ierr);

      // extract preconditioner context
      PC subPc;
      ierr = KSPGetPC(subKsp, &subPc); CHKERRV(ierr);

      // set type of preconditioner
      ierr = PCSetType(subPc, subPcType); CHKERRV(ierr);

      // set Hypre Options from Python config
      if (subPcType == std::string(PCHYPRE))
      {
        std::string hypreOptions = this->specificSettings_.getOptionString("hypreOptions", "");
        PetscOptionsInsertString(NULL, hypreOptions.c_str());
        
        // if one of the hypre preconditioners is in subPreconditionerType, subPcType was set to HYPRE, now set the chosen preconditioner as -sub_pc_hypre_type
        if (subPreconditionerType == "euclid" || subPreconditionerType == "pilut" || subPreconditionerType == "parasails" 
          || subPreconditionerType == "boomeramg" || subPreconditionerType == "ams" || subPreconditionerType == "ads")
        {
#if defined(PETSC_HAVE_HYPRE)
          ierr = PCHYPRESetType(subPc, subPreconditionerType.c_str()); CHKERRV(ierr);
#else
          LOG(ERROR) << "Petsc is not compiled with HYPRE!";
#endif          
          LOG(DEBUG) << "set sub_pc_hypre_type to " << subPreconditionerType;
        }
      }
      
      // set coordinates
      ierr = PCSetCoordinates(subPc, 3, nNodesLocalBlock, nodePositionCoordinatesForPreconditioner.data()); CHKERRV(ierr);
    }
  }

}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusionMuscle,typename FiniteElementMethodDiffusionFat>
void MultidomainWithFatSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle,FiniteElementMethodDiffusionFat>::
updateSystemMatrix()
{
  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_+std::string("_reassemble"));

  // assemble stiffness and mass matrices again
  this->finiteElementMethodDiffusion_.setStiffnessMatrix();
  this->finiteElementMethodDiffusion_.setMassMatrix();
  
  if (useLumpedMassMatrix_)
  {
    this->finiteElementMethodDiffusion_.setInverseLumpedMassMatrix();
  }
  
  this->finiteElementMethodFat_.setStiffnessMatrix();
  this->finiteElementMethodDiffusionTotal_.setStiffnessMatrix();

  LOG(DEBUG) << "rebuild system matrix";
#ifdef DUMP_REBUILT_SYSTEM_MATRIX
  static int counter = 0;
  std::stringstream s; 
  s << counter;

  PetscUtility::dumpMatrix(s.str()+"finiteElementMethodDiffusion_stiffness", "matlab", this->finiteElementMethodDiffusion_.data().stiffnessMatrix()->valuesGlobal(), MPI_COMM_WORLD);
  PetscUtility::dumpMatrix(s.str()+"finiteElementMethodDiffusion_mass", "matlab", this->finiteElementMethodDiffusion_.data().massMatrix()->valuesGlobal(), MPI_COMM_WORLD);
  PetscUtility::dumpMatrix(s.str()+"finiteElementMethodFat_stiffness", "matlab", this->finiteElementMethodFat_.data().stiffnessMatrix()->valuesGlobal(), MPI_COMM_WORLD);
#endif

  // compute new entries for submatrices, except B,C,D and E
  setSystemMatrixSubmatrices(this->timeStepWidthOfSystemMatrix_);
  
  // also compute new entries for the matrices B, C, D and E
  updateBorderMatrices();

  // create the system matrix again
  this->createSystemMatrixFromSubmatrices();

#ifdef DUMP_REBUILT_SYSTEM_MATRIX

  for (int i = 0; i < this->submatricesSystemMatrix_.size(); i++)
  {
    if (this->submatricesSystemMatrix_[i])
    {
      std::stringstream name;
      name << s.str() << "_submatrix_" << i;
      PetscUtility::dumpMatrix(name.str(), "matlab", this->submatricesSystemMatrix_[i], MPI_COMM_WORLD);
    }
  }

  PetscUtility::dumpMatrix(s.str()+"new_system_matrix", "matlab", this->singleSystemMatrix_, MPI_COMM_WORLD);
  counter++;
#endif

  // stop duration measurement
  if (this->durationLogKey_ != "")
  {
    Control::PerformanceMeasurement::stop(this->durationLogKey_+std::string("_reassemble"));

    LOG(INFO) << "Rebuilt multidomain system matrix in " << Control::PerformanceMeasurement::getDuration(this->durationLogKey_+std::string("_reassemble"), false) << "s.";
  }
}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusionMuscle,typename FiniteElementMethodDiffusionFat>
void MultidomainWithFatSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle,FiniteElementMethodDiffusionFat>::
solveLinearSystem()
{
  VLOG(1) << "in solveLinearSystem";

  // configure that the initial value for the iterative solver is the value in solution, not zero
  PetscErrorCode ierr;
  if (this->initialGuessNonzero_)
  {
    LOG(DEBUG) << "set initial guess nonzero";
    ierr = KSPSetInitialGuessNonzero(*this->linearSolver_->ksp(), PETSC_TRUE); CHKERRV(ierr);
  }

  // transform input phi_b to entry in solution vector without shared dofs, this->subvectorsSolution_[this->nCompartments_+1]
  copyPhiBToSolution();

  // compute b_k, the top right hand side entry
  // b_k = b1_ * Vm^(i) + b2_ * phi_e^(i)
  for (int k = 0; k < this->nCompartments_; k++)
  {
    Vec vm_k = this->dataMultidomain_.transmembranePotential(k)->valuesGlobal();     // this is for V_mk^(i)
    Vec phie_k = this->dataMultidomain_.extraCellularPotential()->valuesGlobal();    // this is phi_ek^(i)

    // compute temporary_ = b1_[k]*vm_k
    ierr = MatMult(b1_[k], vm_k, temporary_); CHKERRV(ierr);   // y = Ax

    // compute b_k = temporary_ + b2_[k]*phie_k = b1_[k]*vm_k + b2_[k]*phie_k
    ierr = MatMultAdd(b2_[k], phie_k, temporary_, this->subvectorsRightHandSide_[k]); CHKERRV(ierr);   // v3 = v2 + A * v1, MatMultAdd(Mat mat,Vec v1,Vec v2,Vec v3)
  }

  // copy the values from the nested Petsc Vec nestedRightHandSide_ to the single Vec, singleRightHandSide_, that contains all entries
  NestedMatVecUtility::createVecFromNestedVec(this->nestedRightHandSide_, this->singleRightHandSide_, data().functionSpace()->meshPartition()->rankSubset());

  if (VLOG_IS_ON(1))
  {
    VLOG(1) << "this->nestedRightHandSide_: " << PetscUtility::getStringVector(this->nestedRightHandSide_);
    VLOG(1) << "this->singleRightHandSide_: " << PetscUtility::getStringVector(this->singleRightHandSide_);
  }
  
  bool hasSolverConverged = false;

  // try up to three times to solve the system
  for (int solveNo = 0; solveNo < 5; solveNo++)
  {
    if (solveNo == 0 || solveNo == 1)
    {
      // copy the values from the nested Petsc Vec,nestedSolution_, to the single Vec, singleSolution_, that contains all entries
      NestedMatVecUtility::createVecFromNestedVec(this->nestedSolution_, this->singleSolution_, data().functionSpace()->meshPartition()->rankSubset());
    }

    // Solve the linear system
    // using single Vecs and Mats that contain all values directly  (singleSolution_, singleRightHandSide_, singleSystemMatrix_)
    // This is better compared to using the nested Vec's, because more solvers are available for normal Vec's.
    if (this->showLinearSolverOutput_)
    {
      // solve and show information on convergence
      if (solveNo == 0 || !this->alternativeLinearSolver_)
        hasSolverConverged = this->linearSolver_->solve(this->singleRightHandSide_, this->singleSolution_, "Linear system of multidomain problem solved");
      else
        hasSolverConverged = this->alternativeLinearSolver_->solve(this->singleRightHandSide_, this->singleSolution_, "Linear system of multidomain problem solved");
    }
    else
    {
      // solve without showing output
      if (solveNo == 0 || !this->alternativeLinearSolver_)
        hasSolverConverged = this->linearSolver_->solve(this->singleRightHandSide_, this->singleSolution_);
      else
        hasSolverConverged = this->alternativeLinearSolver_->solve(this->singleRightHandSide_, this->singleSolution_);
    }
    if (hasSolverConverged)
    {
      break;
    }
    else
    {
      if (this->alternativeLinearSolver_)
        LOG(WARNING) << "Solver has not converged, try again with alternative linear solver " << solveNo << "/5";
      else
        LOG(WARNING) << "Solver has not converged, try again " << solveNo << "/5";
    }
  }

  // store the last number of iterations
  this->lastNumberOfIterations_ = this->linearSolver_->lastNumberOfIterations();

  // copy the values back from the single Vec, singleSolution_, that contains all entries 
  // to the nested Petsc Vec, nestedSolution_ which contains the components in subvectorsSolution_
  NestedMatVecUtility::fillNestedVec(this->singleSolution_, this->nestedSolution_);

  // the vector for phi_b in nestedSolution_ contains only entries for non-border dofs, 
  // copy all values and the border dof values to the proper phi_b which is dataFat_.extraCellularPotentialFat()->valuesGlobal()
  copySolutionToPhiB();

  LOG(DEBUG) << "after linear solver:";
  LOG(DEBUG) << "extracellularPotentialFat: " << PetscUtility::getStringVector(dataFat_.extraCellularPotentialFat()->valuesGlobal());
  LOG(DEBUG) << *dataFat_.extraCellularPotentialFat();

}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusionMuscle,typename FiniteElementMethodDiffusionFat>
void MultidomainWithFatSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle,FiniteElementMethodDiffusionFat>::
callOutputWriter(int timeStepNo, double currentTime, int callCountIncrement)
{
  // write current output values
  this->outputWriterManager_.writeOutput(this->dataFat_, timeStepNo, currentTime, callCountIncrement);
}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusionMuscle,typename FiniteElementMethodDiffusionFat>
typename MultidomainWithFatSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle,FiniteElementMethodDiffusionFat>::DataFat &MultidomainWithFatSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle,FiniteElementMethodDiffusionFat>::
data()
{
  return dataFat_;
}

} // namespace TimeSteppingScheme
