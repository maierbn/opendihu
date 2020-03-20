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
  MultidomainSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle>::initialize();

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

  // θ value for Crank-Nicolson scheme
  theta_ = this->specificSettings_.getOptionDouble("theta", 0.5, PythonUtility::Positive);

  DihuContext::solverStructureVisualizer()->setOutputConnectorData(this->getOutputConnectorData());

  // initialize rhs and solution vector, the entries up to this->nCompartments_ were already set by the MultidomainSolver
  this->subvectorsRightHandSide_.resize(this->nCompartments_+2);
  this->subvectorsSolution_.resize(this->nCompartments_+2);

  // resize the submatrices vector
  setSystemMatrixSubmatrices(this->timeStepWidthOfSystemMatrix_);

  // initialize sharedNodes_ 
  findSharedNodesBetweenMuscleAndFat();

  // initialize I_ΓM and -I_ΓM and set in submatrices for system matrix, also set the last zero entry of the rhs
  initializeBorderVariables();

  // create the actual system matrix from the initialize submatrices
  this->createSystemMatrix();

  // system to be solved:
  //
  // [A^1_Vm,Vm   |            |             | B^1_Vm,phie |             ]   [ V^1_m^(i+1) ]    [b^1^(i)]
  // [            | A^2_Vm,Vm  |             | B^2_Vm,phie |             ]   [ V^2_m^(i+1) ]    [b^2^(i)]
  // [   ...      |            | A^M_Vm,Vm   | B^M_Vm,phie |             ] * [ V^M_m^(i+1) ] =  [b^M^(i)]
  // [B^1_phie,Vm |B^2_phie,Vm | B^M_phie,Vm | B_phie,phie | C_phib,phib ]   [ phi_e^(i+1) ]    [ 0     ]
  // [            |            |             | I_ΓM        |  -I_ΓM      ]   [ phi_b^(i+1) ]    [ 0     ]  <- size of last row is number of shared nodes between muscle and fat

  PetscErrorCode ierr;

  // create temporary vector which is needed in computation of rhs
  ierr = MatCreateVecs(b1_[0], &temporary_, NULL); CHKERRV(ierr);

  // initialize subvectors for rhs
  // entry nCompartments was already set to zero by MultidomainSolver, entry nCompartments+1 was already set by initializeBorderVariables()
  for (int k = 0; k < this->nCompartments_; k++)
  {    
    // initialize top part of rhs
    ierr = MatCreateVecs(b1_[k], &this->subvectorsRightHandSide_[k], NULL); CHKERRV(ierr);
  }

  // initialize subvectors for solution
  // entries for Vm from 0 to nCompartments-1 and for phi_e at nCompartments were already set by MultidomainSolver
  // set subvector for phi_b, the body potential
  this->subvectorsSolution_[this->nCompartments_+1]      = dataFat_.extraCellularPotentialFat()->valuesGlobal();    // this is for phi_b
   
  // clear values in solution
  ierr = VecZeroEntries(this->subvectorsSolution_[this->nCompartments_+1]); CHKERRV(ierr);

  // re-create the nested vectors
  LOG(DEBUG) << "create nested vector";
  ierr = VecCreateNest(this->rankSubset_->mpiCommunicator(), this->nCompartments_+2, NULL,
                       this->subvectorsRightHandSide_.data(), &this->nestedRightHandSide_); CHKERRV(ierr);

  ierr = VecCreateNest(this->rankSubset_->mpiCommunicator(), this->nCompartments_+2, NULL, 
                       this->subvectorsSolution_.data(), &this->nestedSolution_); CHKERRV(ierr);

  // copy the values from a nested Petsc Vec to a single Vec that contains all entries
  NestedMatVecUtility::createVecFromNestedVec(this->nestedRightHandSide_, this->singleRightHandSide_, data().functionSpace()->meshPartition()->rankSubset());

}

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

  std::vector<std::vector<Vec3>> nodePositions(2);
  std::vector<std::vector<std::pair<Vec3,node_no_t>>> nodePositionsNodes(2);  // the node positions with nodes for every submesh

  // get geometry fields for function spaces
  functionSpaceMuscle->geometryField().getValuesWithGhosts(nodePositions[0]);
  functionSpaceFat->geometryField().getValuesWithGhosts(nodePositions[1]);


  // iterate over submeshes and save all node positions
  for(int i = 0; i < 2; i++)
  {
    // store node positions with local node nos
    for (node_no_t nodeNoLocal = 0; nodeNoLocal < functionSpaceMuscle->nNodesLocalWithGhosts(); nodeNoLocal++)
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

        break;
      }
    }
  }
}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusionMuscle,typename FiniteElementMethodDiffusionFat>
void MultidomainWithFatSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle,FiniteElementMethodDiffusionFat>::
initializeBorderVariables()
{
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
  int nEntriesLocal = sharedNodes_.size(); 
  ierr = VecSetSizes(vecZero, nEntriesLocal, PETSC_DECIDE); CHKERRV(ierr);

  // set sparsity type and other options
  ierr = VecSetFromOptions(vecZero); CHKERRV(ierr);

  // set all entries to 0.0
  ierr = VecZeroEntries(vecZero); CHKERRV(ierr);

  // store the vector at the bottom of the rhs
  this->subvectorsRightHandSide_[this->nCompartments_+1] = vecZero;

  // -------
  // create matrix I_ΓM (gamma0) with ones for nodes on Γ_M within the muscle domain
  Mat gamma0;
  ierr = MatCreate(mpiCommunicator, &gamma0); CHKERRV(ierr);

  // initialize size
  int nRowsLocal = sharedNodes_.size(); 
  int nColumnsLocalGamma0 = this->finiteElementMethodDiffusion_.functionSpace()->nDofsLocalWithoutGhosts();
  ierr = MatSetSizes(gamma0, nRowsLocal, nColumnsLocalGamma0, PETSC_DECIDE, PETSC_DECIDE); CHKERRV(ierr);
  ierr = MatSetType(gamma0, MATAIJ); CHKERRV(ierr);
  
  // sparse matrix: preallocation of internal data structure
  // http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MATAIJ.html#MATAIJ
  // MATAIJ = "aij" - A matrix type to be used for sparse matrices. This matrix type is identical to MATSEQAIJ when constructed with a single process communicator, and MATMPIAIJ otherwise.
  // As a result, for single process communicators, MatSeqAIJSetPreallocation is supported, and similarly MatMPIAIJSetPreallocation is supported for communicators controlling multiple processes.
  // It is recommended that you call both of the above preallocation routines for simplicity.
  PetscInt diagonalNonZeros = 1;
  PetscInt offdiagonalNonZeros = 1;
  ierr = MatSeqAIJSetPreallocation(gamma0, diagonalNonZeros, NULL); CHKERRV(ierr);
  ierr = MatMPIAIJSetPreallocation(gamma0, diagonalNonZeros, NULL, offdiagonalNonZeros, NULL); CHKERRV(ierr);
  
  // initialize name
  name = "gamma0";
  ierr = PetscObjectSetName((PetscObject)gamma0, name.c_str()); CHKERRV(ierr);

  // -------
  // create matrix I_ΓM (gamma1) with ones for nodes on Γ_M within the fat domain
  Mat gamma1;
  ierr = MatCreate(mpiCommunicator, &gamma1); CHKERRV(ierr);
  
  // initialize size
  int nColumnsLocalGamma1 = finiteElementMethodFat_.functionSpace()->nDofsLocalWithoutGhosts();
  ierr = MatSetSizes(gamma1, nRowsLocal, nColumnsLocalGamma1, PETSC_DECIDE, PETSC_DECIDE); CHKERRV(ierr);
  ierr = MatSetType(gamma1, MATAIJ); CHKERRV(ierr);
  
  // sparse matrix: preallocation of internal data structure
  // http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MATAIJ.html#MATAIJ
  // MATAIJ = "aij" - A matrix type to be used for sparse matrices. This matrix type is identical to MATSEQAIJ when constructed with a single process communicator, and MATMPIAIJ otherwise.
  // As a result, for single process communicators, MatSeqAIJSetPreallocation is supported, and similarly MatMPIAIJSetPreallocation is supported for communicators controlling multiple processes.
  // It is recommended that you call both of the above preallocation routines for simplicity.
  ierr = MatSeqAIJSetPreallocation(gamma1, diagonalNonZeros, NULL); CHKERRV(ierr);
  ierr = MatMPIAIJSetPreallocation(gamma1, diagonalNonZeros, NULL, offdiagonalNonZeros, NULL); CHKERRV(ierr);
  
  // initialize name
  name = "gamma1";
  ierr = PetscObjectSetName((PetscObject)gamma1, name.c_str()); CHKERRV(ierr);

  // set entries of matrices gamma0 and gamma1
  const int nDofsPerNode = this->finiteElementMethodDiffusion_.functionSpace()->nDofsPerNode();
  std::vector<dof_no_t> dofNosMuscle;
  std::vector<dof_no_t> dofNosFat;
  dofNosMuscle.reserve(sharedNodes_.size());
  dofNosFat.reserve(dofNosFat.size());

  // loop over entries <nodeNoMuscle,nodeNoFat> that are shared nodes between the two meshes
  for (std::pair<node_no_t,node_no_t> nodes: sharedNodes_)
  {
    for (int nodalDofNo = 0; nodalDofNo < nDofsPerNode; nodalDofNo++)
    {
      dof_no_t dofNoMuscle = nodes.first*nDofsPerNode + nodalDofNo;
      dof_no_t dofNoFat = nodes.second*nDofsPerNode + nodalDofNo;
      
      dofNosMuscle.push_back(dofNoMuscle);
      dofNosFat.push_back(dofNoFat);
    }
  }

  assert (nRowsLocal == nEntriesLocal);
  std::vector<double> ones(nEntriesLocal, 1.0);
  std::vector<PetscInt> rowIndices(nRowsLocal);
  std::iota(rowIndices.begin(), rowIndices.end(), 0);
  
  ierr = MatSetValues(gamma0, nRowsLocal, rowIndices.data(), nColumnsLocalGamma0, dofNosMuscle.data(), ones.data(), INSERT_VALUES); CHKERRV(ierr);
  ierr = MatSetValues(gamma1, nRowsLocal, rowIndices.data(), nColumnsLocalGamma1, dofNosFat.data(), ones.data(), INSERT_VALUES); CHKERRV(ierr);

  ierr = MatScale(gamma1, -1.0); CHKERRV(ierr);

  // store the matrices in the system matrix
  assert (this->nColumnSubmatricesSystemMatrix_ == this->nCompartments_+2);
  this->submatricesSystemMatrix_[(this->nCompartments_+1)*this->nColumnSubmatricesSystemMatrix_ + this->nCompartments_] = gamma0;
  this->submatricesSystemMatrix_[(this->nCompartments_+1)*this->nColumnSubmatricesSystemMatrix_ + this->nCompartments_ + 1] = gamma1;
}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusionMuscle,typename FiniteElementMethodDiffusionFat>
void MultidomainWithFatSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle,FiniteElementMethodDiffusionFat>::
callOutputWriter(int timeStepNo, double currentTime)
{
  // write current output values
  this->outputWriterManager_.writeOutput(this->dataFat_, timeStepNo, currentTime);
}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusionMuscle,typename FiniteElementMethodDiffusionFat>
void MultidomainWithFatSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle,FiniteElementMethodDiffusionFat>::
setSystemMatrixSubmatrices(double timeStepWidth)
{
  // initialize number of submatrix rows in the system matrix
  this->nColumnSubmatricesSystemMatrix_ = this->nCompartments_+2;

  this->submatricesSystemMatrix_.resize(MathUtility::sqr(this->nColumnSubmatricesSystemMatrix_),NULL);

  // system to be solved:
  //
  // [A^1_Vm,Vm   |            |             | B^1_Vm,phie |             ]   [ V^1_m^(i+1) ]    [b^1^(i)]
  // [            | A^2_Vm,Vm  |             | B^2_Vm,phie |             ]   [ V^2_m^(i+1) ]    [b^2^(i)]
  // [   ...      |            | A^M_Vm,Vm   | B^M_Vm,phie |             ] * [ V^M_m^(i+1) ] =  [b^M^(i)]
  // [B^1_phie,Vm |B^2_phie,Vm | B^M_phie,Vm | B_phie,phie | C_phib,phib ]   [ phi_e^(i+1) ]    [ 0     ]
  // [            |            |             | I_ΓM        |  -I_ΓM      ]   [ phi_b^(i+1) ]    [ 0     ]  <- size of last row is number of shared nodes between muscle and fat

  LOG(TRACE) << "setSystemMatrix";

  // fill this->submatricesSystemMatrix_, empty this->submatricesSystemMatrix_ may be NULL
  // stiffnessMatrix and inverse lumped mass matrix without prefactor
  Mat stiffnessMatrix = this->finiteElementMethodDiffusion_.data().stiffnessMatrix()->valuesGlobal();
  Mat massMatrix = this->finiteElementMethodDiffusion_.data().massMatrix()->valuesGlobal();


  PetscErrorCode ierr;
  // set all this->submatricesSystemMatrix_
  for (int k = 0; k < this->nCompartments_; k++)
  {
    // right column matrix
    double prefactor = theta_ / (this->am_[k]*this->cm_[k]);

    VLOG(2) << "k=" << k << ", am: " << this->am_[k] << ", cm: " << this->cm_[k] << ", prefactor: " << prefactor;

    // create matrix as theta/(Am*Cm)*K
    Mat matrixOnRightColumn;
    ierr = MatConvert(stiffnessMatrix, MATSAME, MAT_INITIAL_MATRIX, &matrixOnRightColumn); CHKERRV(ierr);

    // scale matrix on right column with prefactor
    ierr = MatScale(matrixOnRightColumn, prefactor); CHKERRV(ierr);

    // set on right column of the system matrix
    this->submatricesSystemMatrix_[k*this->nColumnSubmatricesSystemMatrix_ + this->nCompartments_] = matrixOnRightColumn;

    // ---
    // diagonal matrix
    // copy right block matrix also to diagonal matrix
    Mat matrixOnDiagonalBlock;
    ierr = MatConvert(matrixOnRightColumn, MATSAME, MAT_INITIAL_MATRIX, &matrixOnDiagonalBlock); CHKERRV(ierr);

    // add scaled mass matrix, -1/dt*M,  AXPY: Y = a*X + Y  MatAXPY(Y,a,X,SAME_NONZERO_PATTERN)
    prefactor = -1/timeStepWidth;
    ierr = MatAXPY(matrixOnDiagonalBlock, prefactor, massMatrix, SAME_NONZERO_PATTERN); CHKERRV(ierr);

    // set on diagonal
    this->submatricesSystemMatrix_[k*this->nColumnSubmatricesSystemMatrix_ + k] = matrixOnDiagonalBlock;

    // ---
    // bottom row matrices
    // stiffnessMatrixWithPrefactor is f_k*K
    Mat stiffnessMatrixWithPrefactor = this->finiteElementMethodDiffusionCompartment_[k].data().stiffnessMatrix()->valuesGlobal();

    // create matrix as copy of stiffnessMatrix
    Mat matrixOnBottomRow;
    ierr = MatConvert(stiffnessMatrixWithPrefactor, MATSAME, MAT_INITIAL_MATRIX, &matrixOnBottomRow); CHKERRV(ierr);

    // set on bottom row of the system matrix
    this->submatricesSystemMatrix_[this->nCompartments_*this->nColumnSubmatricesSystemMatrix_ + k] = matrixOnBottomRow;
  }

  // set bottom right matrix 1
  Mat stiffnessMatrixBottomRightMuscle = this->finiteElementMethodDiffusionTotal_.data().stiffnessMatrix()->valuesGlobal();

  // set on bottom right
  this->submatricesSystemMatrix_[this->nCompartments_*this->nColumnSubmatricesSystemMatrix_ + this->nCompartments_] = stiffnessMatrixBottomRightMuscle;
  
  // set bottom right matrix 2
  Mat stiffnessMatrixBottomRightFat = this->finiteElementMethodFat_.data().stiffnessMatrix()->valuesGlobal();

  // set on bottom right
  this->submatricesSystemMatrix_[this->nCompartments_*this->nColumnSubmatricesSystemMatrix_ + this->nCompartments_+1] = stiffnessMatrixBottomRightFat;

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

    // add scaled mass matrix, -1/dt*M,  AXPY: Y = a*X + Y  MatAXPY(Y,a,X,SAME_NONZERO_PATTERN)
    prefactor = -1/timeStepWidth;
    ierr = MatAXPY(b1_[k], prefactor, massMatrix, SAME_NONZERO_PATTERN); CHKERRV(ierr);


    // set b2_ = (θ-1)*K_sigmai^k
    ierr = MatConvert(stiffnessMatrix, MATSAME, MAT_INITIAL_MATRIX, &b2_[k]); CHKERRV(ierr);
    
    prefactor = theta_ - 1;
    ierr = MatScale(b2_[k], prefactor); CHKERRV(ierr);
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

  // copy the values from a nested Petsc Vec to a single Vec that contains all entries
  NestedMatVecUtility::createVecFromNestedVec(this->nestedRightHandSide_, this->singleRightHandSide_, data().functionSpace()->meshPartition()->rankSubset());

  // copy the values from a nested Petsc Vec to a single Vec that contains all entries
  NestedMatVecUtility::createVecFromNestedVec(this->nestedSolution_, this->singleSolution_, data().functionSpace()->meshPartition()->rankSubset());

  // solve the linear system
  // this can be done using the nested Vecs and nested Mat (nestedSolution_, nestedRightHandSide_, nestedSystemMatrix_),
  // or the single Vecs and Mats that contain all values directly  (singleSolution_, singleRightHandSide_, singleSystemMatrix_) 


  if (this->showLinearSolverOutput_)
  {
    this->linearSolver_->solve(this->singleRightHandSide_, this->singleSolution_, "Linear system of multidomain problem solved");
  }
  else
  {
    // solve without showing output
    this->linearSolver_->solve(this->singleRightHandSide_, this->singleSolution_);
  }
  
  // copy the values back from a single Vec that contains all entries to a nested Petsc Vec
  NestedMatVecUtility::fillNestedVec(this->singleSolution_, this->nestedSolution_);
}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusionMuscle,typename FiniteElementMethodDiffusionFat>
typename MultidomainWithFatSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle,FiniteElementMethodDiffusionFat>::DataFat &MultidomainWithFatSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle,FiniteElementMethodDiffusionFat>::
data()
{
  return dataFat_;
}

} // namespace TimeSteppingScheme
