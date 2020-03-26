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

  DihuContext::solverStructureVisualizer()->setOutputConnectorData(this->getOutputConnectorData());

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

  // set matrix used for linear solver and preconditioner to ksp context
  assert(this->linearSolver_->ksp());
  PetscErrorCode ierr;
  ierr = KSPSetOperators(*this->linearSolver_->ksp(), this->singleSystemMatrix_, this->singleSystemMatrix_); CHKERRV(ierr);

  // set the nullspace of the matrix
  // as we have Neumann boundary conditions, constant functions are in the nullspace of the matrix
  MatNullSpace nullSpace;
  ierr = MatNullSpaceCreate(data().functionSpace()->meshPartition()->mpiCommunicator(), PETSC_TRUE, 0, PETSC_NULL, &nullSpace); CHKERRV(ierr);
  ierr = MatSetNullSpace(this->singleSystemMatrix_, nullSpace); CHKERRV(ierr);
  ierr = MatSetNearNullSpace(this->singleSystemMatrix_, nullSpace); CHKERRV(ierr); // for multigrid methods
  //ierr = MatNullSpaceDestroy(&nullSpace); CHKERRV(ierr);

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

    // matrix B
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

  // copy the values from the nested Petsc Vec,nestedSolution_, to the single Vec, singleSolution_, that contains all entries
  NestedMatVecUtility::createVecFromNestedVec(this->nestedSolution_, this->singleSolution_, data().functionSpace()->meshPartition()->rankSubset());

  if (VLOG_IS_ON(1))
  {
    VLOG(1) << "this->nestedRightHandSide_: " << PetscUtility::getStringVector(this->nestedRightHandSide_);
    VLOG(1) << "this->singleRightHandSide_: " << PetscUtility::getStringVector(this->singleRightHandSide_);
  }

  // Solve the linear system
  // using single Vecs and Mats that contain all values directly  (singleSolution_, singleRightHandSide_, singleSystemMatrix_)
  // This is better compared to using the nested Vec's, because more solvers are available for normal Vec's.
  if (this->showLinearSolverOutput_)
  {
    // solve and show information on convergence
    this->linearSolver_->solve(this->singleRightHandSide_, this->singleSolution_, "Linear system of multidomain problem solved");
  }
  else
  {
    // solve without showing output
    this->linearSolver_->solve(this->singleRightHandSide_, this->singleSolution_);
  }
  
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
callOutputWriter(int timeStepNo, double currentTime)
{
  // write current output values
  this->outputWriterManager_.writeOutput(this->dataFat_, timeStepNo, currentTime);
}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusionMuscle,typename FiniteElementMethodDiffusionFat>
typename MultidomainWithFatSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle,FiniteElementMethodDiffusionFat>::DataFat &MultidomainWithFatSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle,FiniteElementMethodDiffusionFat>::
data()
{
  return dataFat_;
}

} // namespace TimeSteppingScheme
