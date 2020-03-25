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
  // indicate in solverStructureVisualizer that now a child solver will be initialized
  DihuContext::solverStructureVisualizer()->beginChild("Fat");

  LOG(DEBUG) << "initialize fat FEM";
  // initialize the potential flow finite element method, this also creates the function space
  finiteElementMethodFat_.initialize();

  // indicate in solverStructureVisualizer that the child solver initialization is done
  DihuContext::solverStructureVisualizer()->endChild();

  MultidomainSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle>::initialize();

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
  // [B^1_phie,Vm |B^2_phie,Vm | B^M_phie,Vm | B_phie,phie |      D      ]   [ phi_e^(i+1) ]    [ 0     ]
  // [            |            |             |     E       | C_phib,phib ]   [ phi_b^(i+1) ]    [ 0     ]

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
callOutputWriter(int timeStepNo, double currentTime)
{
  // write current output values
  this->outputWriterManager_.writeOutput(this->dataFat_, timeStepNo, currentTime);
}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusionMuscle,typename FiniteElementMethodDiffusionFat>
void MultidomainWithFatSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle,FiniteElementMethodDiffusionFat>::
setSystemMatrixSubmatrices(double timeStepWidth)
{
  assert(this->finiteElementMethodDiffusionTotal_.data().stiffnessMatrix());
  assert(this->finiteElementMethodDiffusion_.data().stiffnessMatrix());
  assert(this->finiteElementMethodFat_.data().stiffnessMatrix());

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

  PetscErrorCode ierr;
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

    // set on right column of the system matrix
    this->submatricesSystemMatrix_[k*this->nColumnSubmatricesSystemMatrix_ + this->nCompartments_] = matrixOnRightColumn;

    // ---
    // matrix on diagonal, A
    // copy right block matrix also to diagonal matrix
    Mat matrixOnDiagonalBlock;
    ierr = MatConvert(matrixOnRightColumn, MATSAME, MAT_INITIAL_MATRIX, &matrixOnDiagonalBlock); CHKERRV(ierr);

    // add scaled mass matrix, -1/dt*M,  AXPY: Y = a*X + Y  MatAXPY(Y,a,X,SAME_NONZERO_PATTERN)
    prefactor = -1/timeStepWidth;
    ierr = MatAXPY(matrixOnDiagonalBlock, prefactor, massMatrix, SAME_NONZERO_PATTERN); CHKERRV(ierr);

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

  // set bottom right matrix, B
  Mat stiffnessMatrixBottomRightMuscle = this->finiteElementMethodDiffusionTotal_.data().stiffnessMatrix()->valuesGlobal();

  // set on bottom right
  this->submatricesSystemMatrix_[this->nCompartments_*this->nColumnSubmatricesSystemMatrix_ + this->nCompartments_] = stiffnessMatrixBottomRightMuscle;
  
  // set bottom right matrix, C
  Mat stiffnessMatrixBottomRightFat = this->finiteElementMethodFat_.data().stiffnessMatrix()->valuesGlobal();

  // set on bottom right
  this->submatricesSystemMatrix_[(this->nCompartments_+1)*this->nColumnSubmatricesSystemMatrix_ + this->nCompartments_+1] = stiffnessMatrixBottomRightFat;

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
