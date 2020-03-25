#include "specialized_solver/multidomain_solver/multidomain_with_fat_solver.h"

#include <Python.h>  // has to be the first included header

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
        
        // add dofs at node to set of border nodes of muscle mesh
        for (int i = 0; i < FunctionSpace::nDofsPerNode(); i++)
        {
          dof_no_t dofNoLocal = FunctionSpace::nDofsPerNode()*nodeNoLocal + i;
          borderDofsMuscle_.insert(dofNoLocal);
        }

        // add dofs at node to set of border nodes of fat mesh
        for (int i = 0; i < FunctionSpace::nDofsPerNode(); i++)
        {
          dof_no_t dofNoLocal = FunctionSpace::nDofsPerNode()*nodeNoLocalOtherMesh + i;
          borderDofsFat_.insert(dofNoLocal);
        }

        break;
      }
    }
  }

  LOG(DEBUG) << "functionSpaceMuscle: " << *functionSpaceMuscle->meshPartition();
  LOG(DEBUG) << "functionSpaceFat: " << *functionSpaceFat->meshPartition();

  LOG(DEBUG) << "sharedNodes_:\n" << sharedNodes_;
  LOG(DEBUG) << "borderDofsMuscle_:\n" << borderDofsMuscle_;
  LOG(DEBUG) << "borderDofsFat_:\n" << borderDofsFat_;

  // fill variables borderElementDofsMuscle_ and borderElementDofsFat_

  const int nDofsPerElement = FunctionSpace::nDofsPerElement();

  // loop over elements of muscle mesh
  for (element_no_t elementNoLocal = 0; elementNoLocal < functionSpaceMuscle->nElementsLocal(); elementNoLocal++)
  {
    // get indices of element-local dofs
    std::array<dof_no_t,nDofsPerElement> dofNosLocal = functionSpaceMuscle->getElementDofNosLocal(elementNoLocal);

    // check if node indices are part of sharedNodes_
    int dofIndex = 0;
    for (int elementalDofIndex = 0; elementalDofIndex < nDofsPerElement; elementalDofIndex++)
    {
      dof_no_t dofNoLocal = dofNosLocal[elementalDofIndex];
      if (borderDofsMuscle_.find(dofNoLocal) != borderDofsMuscle_.end())
      {
        if (borderElementDofsMuscle_.find(elementNoLocal) == borderElementDofsMuscle_.end())
        {
          borderElementDofsMuscle_[elementNoLocal] = std::array<dof_no_t,::FunctionSpace::FunctionSpaceBaseDim<2,typename FunctionSpace::BasisFunction>::nDofsPerElement()>();
        }
        LOG(DEBUG) << "set border element " << elementNoLocal << ", dofIndex " << dofIndex << " to " << elementalDofIndex;
        borderElementDofsMuscle_[elementNoLocal][dofIndex++] = elementalDofIndex;
      }
    }
  }
  // loop over elements of fat mesh
  for (element_no_t elementNoLocal = 0; elementNoLocal < functionSpaceFat->nElementsLocal(); elementNoLocal++)
  {
    // get indices of element-local dofs
    std::array<dof_no_t,nDofsPerElement> dofNosLocal = functionSpaceFat->getElementDofNosLocal(elementNoLocal);

    // check if node indices are part of sharedNodes_
    int dofIndex = 0;
    for (int elementalDofIndex = 0; elementalDofIndex < nDofsPerElement; elementalDofIndex++)
    {
      dof_no_t dofNoLocal = dofNosLocal[elementalDofIndex];
      if (borderDofsFat_.find(dofNoLocal) != borderDofsFat_.end())
      {
        borderElementDofsFat_[elementNoLocal][dofIndex++] = elementalDofIndex;
      }
    }
  }
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
  PetscInt nNonZerosDiagonal = 1;
  PetscInt nNonZerosOffdiagonal = 1;
  ierr = MatSeqAIJSetPreallocation(gamma0, nNonZerosDiagonal, NULL); CHKERRV(ierr);
  ierr = MatMPIAIJSetPreallocation(gamma0, nNonZerosDiagonal, NULL, nNonZerosOffdiagonal, NULL); CHKERRV(ierr);
  
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
  ierr = MatSeqAIJSetPreallocation(gamma1, nNonZerosDiagonal, NULL); CHKERRV(ierr);
  ierr = MatMPIAIJSetPreallocation(gamma1, nNonZerosDiagonal, NULL, nNonZerosOffdiagonal, NULL); CHKERRV(ierr);
  
  // initialize name
  name = "gamma1";
  ierr = PetscObjectSetName((PetscObject)gamma1, name.c_str()); CHKERRV(ierr);

  // set entries of matrices gamma0 and gamma1
  const int nDofsPerNode = this->finiteElementMethodDiffusion_.functionSpace()->nDofsPerNode();

  PetscInt rowNoLocal = 0;
  ierr = MatGetOwnershipRange(gamma0, &rowNoLocal, NULL); CHKERRV(ierr);
  
  // loop over entries <nodeNoMuscle,nodeNoFat> that are shared nodes between the two meshes
  for (std::pair<node_no_t,node_no_t> nodes: sharedNodes_)
  {
    for (int nodalDofNo = 0; nodalDofNo < nDofsPerNode; nodalDofNo++)
    {
      dof_no_t dofNoMuscle = nodes.first*nDofsPerNode + nodalDofNo;
      dof_no_t dofNoFat = nodes.second*nDofsPerNode + nodalDofNo;
      
      ierr = MatSetValue(gamma0, rowNoLocal, dofNoMuscle, 1.0, INSERT_VALUES); CHKERRV(ierr);
      ierr = MatSetValue(gamma1, rowNoLocal, dofNoFat, -1.0, INSERT_VALUES); CHKERRV(ierr);

      rowNoLocal++;
    }
  }


  PetscInt nRowsGamma0, nColumnsGamma0;
  ierr = MatGetSize(gamma0, &nRowsGamma0, &nColumnsGamma0); CHKERRV(ierr);
  PetscInt nRowsGamma1, nColumnsGamma1;
  ierr = MatGetSize(gamma1, &nRowsGamma1, &nColumnsGamma1); CHKERRV(ierr);
  LOG(DEBUG) << "create gamma matrices: gamma0: " << nRowsGamma0 << "x" << nColumnsGamma0 << ", gamma1: " << nRowsGamma1 << "x" << nColumnsGamma1
    << " set at (" << (this->nCompartments_+2) << "," << this->nCompartments_ << ") and (" << (this->nCompartments_+2) << "," << this->nCompartments_+1 << ")";

  // store the matrices in the system matrix
  assert (this->nColumnSubmatricesSystemMatrix_ == this->nCompartments_+3);
  this->submatricesSystemMatrix_[(this->nCompartments_+2)*this->nColumnSubmatricesSystemMatrix_ + this->nCompartments_] = gamma0;
  this->submatricesSystemMatrix_[(this->nCompartments_+2)*this->nColumnSubmatricesSystemMatrix_ + this->nCompartments_ + 1] = gamma1;

  // ---
  // create matrices B_ΓM (bGamma0) and -B_ΓM (bGamma1)
  Mat bGamma0;
  ierr = MatCreate(mpiCommunicator, &bGamma0); CHKERRV(ierr);
  
  // initialize size
  int nRowsLocalBGamma0 = this->finiteElementMethodDiffusion_.functionSpace()->nDofsLocalWithoutGhosts();
  int nColumnsLocalBGamma0 = sharedNodes_.size();
  ierr = MatSetSizes(bGamma0, nRowsLocalBGamma0, nColumnsLocalBGamma0, PETSC_DECIDE, PETSC_DECIDE); CHKERRV(ierr);
  ierr = MatSetType(bGamma0, MATAIJ); CHKERRV(ierr);

  // initialize memory with number of nonzeros
  const int nDofsPerBasis = ::FunctionSpace::FunctionSpaceBaseDim<1,typename FunctionSpace::BasisFunction>::nDofsPerElement();
  const int nOverlaps = (nDofsPerBasis*2 - 1) * nDofsPerNode;   // number of nodes of 2 neighbouring 1D elements (=number of ansatz functions in support of center ansatz function)
  nNonZerosDiagonal = std::pow(nOverlaps, 2);
  nNonZerosOffdiagonal = nNonZerosDiagonal;
    
  ierr = MatSeqAIJSetPreallocation(bGamma0, nNonZerosDiagonal, NULL); CHKERRV(ierr);
  ierr = MatMPIAIJSetPreallocation(bGamma0, nNonZerosDiagonal, NULL, nNonZerosOffdiagonal, NULL); CHKERRV(ierr);
  
  // initialize name
  name = "bGamma0";
  ierr = PetscObjectSetName((PetscObject)bGamma0, name.c_str()); CHKERRV(ierr);

  // create matrix -B_ΓM
  Mat bGamma1;
  ierr = MatCreate(mpiCommunicator, &bGamma1); CHKERRV(ierr);
  
  // initialize size
  int nRowsLocalBGamma1 = this->finiteElementMethodFat_.functionSpace()->nDofsLocalWithoutGhosts();
  ierr = MatSetSizes(bGamma1, nRowsLocalBGamma1, nColumnsLocalBGamma0, PETSC_DECIDE, PETSC_DECIDE); CHKERRV(ierr);
  ierr = MatSetType(bGamma1, MATAIJ); CHKERRV(ierr);
  ierr = MatSeqAIJSetPreallocation(bGamma1, nNonZerosDiagonal, NULL); CHKERRV(ierr);
  ierr = MatMPIAIJSetPreallocation(bGamma1, nNonZerosDiagonal, NULL, nNonZerosOffdiagonal, NULL); CHKERRV(ierr);
  
  // initialize name
  name = "bGamma1";
  ierr = PetscObjectSetName((PetscObject)bGamma1, name.c_str()); CHKERRV(ierr);


  LOG(DEBUG) << "nNonZerosDiagonal: " << nNonZerosDiagonal << ", nNonZerosOffdiagonal: " << nNonZerosOffdiagonal;

  // compute B_ΓM and -B_ΓM (note that bGamma0 != -bGamma1, they use different meshes)
  computeBorderMatrices(bGamma0, bGamma1);

  // store in submatrices of system matrix
  this->submatricesSystemMatrix_[(this->nCompartments_)*this->nColumnSubmatricesSystemMatrix_ + this->nCompartments_ + 2] = bGamma0;
  this->submatricesSystemMatrix_[(this->nCompartments_+1)*this->nColumnSubmatricesSystemMatrix_ + this->nCompartments_ + 2] = bGamma1;
}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusionMuscle,typename FiniteElementMethodDiffusionFat>
void MultidomainWithFatSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusionMuscle,FiniteElementMethodDiffusionFat>::
computeBorderMatrices(Mat bGamma0, Mat bGamma1)
{
  using MuscleFunctionSpace = typename FiniteElementMethodDiffusionMuscle::FunctionSpace;
  using FatFunctionSpace = typename FiniteElementMethodDiffusionFat::FunctionSpace;
  std::shared_ptr<MuscleFunctionSpace> functionSpaceMuscle = this->finiteElementMethodDiffusion_.functionSpace();
  std::shared_ptr<FatFunctionSpace> functionSpaceFat = finiteElementMethodFat_.functionSpace();
  
  // define shortcuts for integrator and basis
  typedef Quadrature::TensorProduct<2,Quadrature::Gauss<3>> Quadrature2D;
  const int nDofsPerElement = FunctionSpace::nDofsPerElement();
  const int nDofsPerElement2D = ::FunctionSpace::FunctionSpaceBaseDim<2,typename FunctionSpace::BasisFunction>::nDofsPerElement();
  typedef MathUtility::Matrix<nDofsPerElement,nDofsPerElement2D> EvaluationsType;
  typedef std::array<
            EvaluationsType,
            Quadrature2D::numberEvaluations()
          > EvaluationsArrayType;     // evaluations[nGP^D][nDofs][nDofs]

  // ---------------------------------------------------------------
  // bGamma0_kj = ∫_ΓM ψj*φk dx, loop ever elements of muscle mesh
  // initialize values to zero
  PetscErrorCode ierr;
  // loop over elements of muscle mesh
  for (element_no_t elementNo = 0; elementNo < functionSpaceMuscle->nElementsLocal(); elementNo++)
  {
    if (borderElementDofsMuscle_.find(elementNo) != borderElementDofsMuscle_.end())
    {
      std::array<dof_no_t,nDofsPerElement> dofNosLocal = functionSpaceMuscle->getElementDofNosLocal(elementNo);

      // loop over all dofs in the element
      for (int i = 0; i < nDofsPerElement; i++)
      {
        // loop over the dofs in the element that are on the border
        for (int jIndex = 0; jIndex < nDofsPerElement2D; jIndex++)
        {
          int j = borderElementDofsMuscle_[elementNo][jIndex];
          dof_no_t dofNoLocalJ = dofNosLocal[j];

          LOG(DEBUG) << "bGamma0(" << dofNosLocal[i] << "," << dofNoLocalJ << ") = 0";
          ierr = MatSetValue(bGamma0, dofNosLocal[i], dofNoLocalJ, 0, INSERT_VALUES); CHKERRV(ierr);
        }
      }
    }
  }

  ierr = MatAssemblyBegin(bGamma0, MAT_FLUSH_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(bGamma0, MAT_FLUSH_ASSEMBLY); CHKERRV(ierr);

  // setup arrays used for integration
  std::array<Vec2, Quadrature2D::numberEvaluations()> samplingPoints = Quadrature2D::samplingPoints();
  EvaluationsArrayType evaluationsArray{};

  // set entries in the matrix
  // loop over elements
  for (element_no_t elementNo = 0; elementNo < functionSpaceMuscle->nElementsLocal(); elementNo++)
  {
    if (borderElementDofsMuscle_.find(elementNo) == borderElementDofsMuscle_.end())
      continue;

    // get the face
    Mesh::face_t face = functionSpaceMuscle->getFaceFromElementalDofNos(borderElementDofsMuscle_[elementNo]);

    LOG(DEBUG) << "In computation of bGamma0, element " << elementNo << ", face " << Mesh::getString(face);

    // get indices of element-local dofs
    std::array<dof_no_t,nDofsPerElement> dofNosLocal = functionSpaceMuscle->getElementDofNosLocal(elementNo);

    // get geometry field (which are the node positions for Lagrange basis and node positions and derivatives for Hermite)
    std::array<Vec3,FunctionSpace::nDofsPerElement()> geometry;
    functionSpaceMuscle->getElementGeometry(elementNo, geometry);

    // compute integral
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPoints.size(); samplingPointIndex++)
    {
      // evaluate function to integrate at samplingPoints[i*2], write value to evaluations[i]
      Vec2 xiSurface = samplingPoints[samplingPointIndex];
      Vec3 xi = Mesh::getXiOnFace(face, xiSurface);

      // compute the 3xD jacobian of the parameter space to world space mapping
      auto jacobian = FunctionSpace::computeJacobian(geometry, xi);
      
      EvaluationsType evaluations;

      // get the factor in the integral that arises from the change in integration domain from world to coordinate space
      double integrationFactor = MathUtility::computeIntegrationFactor<3>(jacobian);

      // loop over pairs of basis functions and evaluation integrand at xi
      // loop over all dofs in the element
      for (int i = 0; i < nDofsPerElement; i++)
      {
        // loop over the dofs in the element that are on the border
        for (int jIndex = 0; jIndex < nDofsPerElement2D; jIndex++)
        {
          int j = borderElementDofsMuscle_[elementNo][jIndex];
          
          double integrand = FunctionSpace::phi(i,xi) * FunctionSpace::phi(j,xi) * integrationFactor;
          
          evaluationsArray[samplingPointIndex](i, j) = integrand;
        }
      }
    }  // function evaluations

    // integrate all values for the (i,j) dof pairs at once
    EvaluationsType integratedValues = Quadrature2D::computeIntegral(evaluationsArray);

    // perform integration and add to entry in rhs vector
    for (int i = 0; i < nDofsPerElement; i++)
    {
      for (int jIndex = 0; jIndex < nDofsPerElement2D; jIndex++)
      {
        int j = borderElementDofsMuscle_[elementNo][jIndex];
        
        // integrate value and set entry in discretization matrix
        double integratedValue = integratedValues(i, j);

        ierr = MatSetValue(bGamma0, dofNosLocal[i], dofNosLocal[j], integratedValue, ADD_VALUES); CHKERRV(ierr);
      }  // j
    }  // i
  }  // elementNo

  // merge local changes in parallel and assemble the matrix (MatAssemblyBegin, MatAssemblyEnd)
  ierr = MatAssemblyBegin(bGamma0, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(bGamma0, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);

  // ---------------------------------------------------------------
  // bGamma1_kj = -∫_ΓM ψj*φk dx, loop ever elements of fat mesh
  // initialize values to zero
  // loop over elements of muscle mesh
  for (element_no_t elementNo = 0; elementNo < functionSpaceFat->nElementsLocal(); elementNo++)
  {
    if (borderElementDofsFat_.find(elementNo) != borderElementDofsFat_.end())
    {
      std::array<dof_no_t,nDofsPerElement> dofNosLocal = functionSpaceFat->getElementDofNosLocal(elementNo);

      // loop over all dofs in the element
      for (int i = 0; i < nDofsPerElement; i++)
      {
        // loop over the dofs in the element that are on the border
        for (int jIndex = 0; jIndex < nDofsPerElement2D; jIndex++)
        {
          int j = borderElementDofsFat_[elementNo][jIndex];
          dof_no_t dofNoLocalJ = dofNosLocal[j];

          ierr = MatSetValue(bGamma1, dofNosLocal[i], dofNoLocalJ, 0, INSERT_VALUES); CHKERRV(ierr);
        }
      }
    }
  }

  ierr = MatAssemblyBegin(bGamma1, MAT_FLUSH_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(bGamma1, MAT_FLUSH_ASSEMBLY); CHKERRV(ierr);

  // set entries in the matrix
  // loop over elements
  for (element_no_t elementNo = 0; elementNo < functionSpaceFat->nElementsLocal(); elementNo++)
  {
    if (borderElementDofsFat_.find(elementNo) == borderElementDofsFat_.end())
      continue;

    // get the face
    Mesh::face_t face = functionSpaceFat->getFaceFromElementalDofNos(borderElementDofsFat_[elementNo]);

    LOG(DEBUG) << "In computation of bGamma1, element " << elementNo << ", face " << Mesh::getString(face);

    // get indices of element-local dofs
    std::array<dof_no_t,nDofsPerElement> dofNosLocal = functionSpaceFat->getElementDofNosLocal(elementNo);

    // get geometry field (which are the node positions for Lagrange basis and node positions and derivatives for Hermite)
    std::array<Vec3,FunctionSpace::nDofsPerElement()> geometry;
    functionSpaceFat->getElementGeometry(elementNo, geometry);

    // compute integral
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPoints.size(); samplingPointIndex++)
    {
      // evaluate function to integrate at samplingPoints[i*2], write value to evaluations[i]
      Vec2 xiSurface = samplingPoints[samplingPointIndex];
      Vec3 xi = Mesh::getXiOnFace(face, xiSurface);

      // compute the 3xD jacobian of the parameter space to world space mapping
      auto jacobian = FunctionSpace::computeJacobian(geometry, xi);
      
      EvaluationsType evaluations;

      // get the factor in the integral that arises from the change in integration domain from world to coordinate space
      double integrationFactor = MathUtility::computeIntegrationFactor<3>(jacobian);

      // loop over pairs of basis functions and evaluation integrand at xi
      // loop over all dofs in the element
      for (int i = 0; i < nDofsPerElement; i++)
      {
        // loop over the dofs in the element that are on the border
        for (int jIndex = 0; jIndex < nDofsPerElement2D; jIndex++)
        {
          int j = borderElementDofsFat_[elementNo][jIndex];
          
          double integrand = FunctionSpace::phi(i,xi) * FunctionSpace::phi(j,xi) * integrationFactor;
          
          evaluationsArray[samplingPointIndex](i, j) = integrand;
        }
      }
    }  // function evaluations

    // integrate all values for the (i,j) dof pairs at once
    EvaluationsType integratedValues = Quadrature2D::computeIntegral(evaluationsArray);

    // perform integration and add to entry in rhs vector
    for (int i = 0; i < nDofsPerElement; i++)
    {
      for (int jIndex = 0; jIndex < nDofsPerElement2D; jIndex++)
      {
        int j = borderElementDofsFat_[elementNo][jIndex];
        
        // integrate value and set entry in discretization matrix
        double integratedValue = -integratedValues(i, j);

        ierr = MatSetValue(bGamma1, dofNosLocal[i], dofNosLocal[j], integratedValue, ADD_VALUES); CHKERRV(ierr);
      }  // j
    }  // i
  }  // elementNo

  // merge local changes in parallel and assemble the matrix (MatAssemblyBegin, MatAssemblyEnd)
  ierr = MatAssemblyBegin(bGamma1, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(bGamma1, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
}


} // namespace TimeSteppingScheme
