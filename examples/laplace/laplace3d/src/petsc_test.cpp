#include <iostream>
#include <cstdlib>

#include "opendihu.h"
#include "partition/partitioned_petsc_vec/partitioned_petsc_vec.h"

int testPlain(int argc, char *argv[])
{

  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  // PetscInitialize(&argc,&argv, NULL, NULL);

  PetscErrorCode ierr;
  //int ownRankNo = 0;
  int ownRankNo = settings.ownRankNo();
  MPI_Comm mpiCommunicator = MPI_COMM_WORLD;
  //MPI_Comm_rank(mpiCommunicator, &ownRankNo);

  int nDofsGlobal = 20;
  int nDofsLocalWithoutGhosts = 0;
  dof_no_t nGhostDofs = 0;
  std::vector<int> localIndices;

  MPI_Barrier(mpiCommunicator);
  LOG(DEBUG) << "ownRankNo: " << ownRankNo;

  LOG(INFO) << "-------------------- START TEST ------------------------";
  Vec rhs;
  /*
   *  3D example, 4 elements, 20 dofs
   *  p0  p1  p2
   * o-o- o- o-o
   * | |  |  | |
   * o-o- o- o-o
   *
   * */

  // create vector
  if (ownRankNo == 0)
  {
    localIndices = std::vector<int>({0,1,2,3,4,5,6,7});
    // create Petsc vecs
    nGhostDofs = 4;
    nDofsLocalWithoutGhosts = 8;

    std::array<int,4> ghostDofNosGlobalPetsc = {8, 9, 10, 11};

    ierr = VecCreateGhost(mpiCommunicator, nDofsLocalWithoutGhosts,
                          nDofsGlobal, nGhostDofs, ghostDofNosGlobalPetsc.data(), &rhs); CHKERRQ(ierr);
                          /*
    ierr = VecCreate(mpiCommunicator, &rhs); CHKERRQ(ierr);
    ierr = VecSetSizes(rhs, nDofsLocalWithoutGhosts, nDofsGlobal); CHKERRQ(ierr);*/
  }
  else if (ownRankNo == 1)
  {
    localIndices = std::vector<int>({8,9,10,11});

    // create Petsc vecs
    nGhostDofs = 4;
    nDofsLocalWithoutGhosts = 4;

    std::array<int,4> ghostDofNosGlobalPetsc = {12, 13, 14, 15};

    ierr = VecCreateGhost(mpiCommunicator, nDofsLocalWithoutGhosts,
                          nDofsGlobal, nGhostDofs, ghostDofNosGlobalPetsc.data(), &rhs); CHKERRQ(ierr);
                          /*
    ierr = VecCreate(mpiCommunicator, &rhs); CHKERRQ(ierr);
    ierr = VecSetSizes(rhs, nDofsLocalWithoutGhosts, nDofsGlobal); CHKERRQ(ierr);*/
  }
  else if (ownRankNo == 2)
  {
    localIndices = std::vector<int>({12,13,14,15,16,17,18,19});

    // create Petsc vecs
    nGhostDofs = 0;
    nDofsLocalWithoutGhosts = 8;

    int ghostDofNosGlobalPetsc = 0;

    ierr = VecCreateGhost(mpiCommunicator, nDofsLocalWithoutGhosts,
                          nDofsGlobal, nGhostDofs, &ghostDofNosGlobalPetsc, &rhs); CHKERRQ(ierr);
                          /*
    ierr = VecCreate(mpiCommunicator, &rhs); CHKERRQ(ierr);
    ierr = VecSetSizes(rhs, nDofsLocalWithoutGhosts, nDofsGlobal); CHKERRQ(ierr);*/
  }

  LOG(DEBUG) << "nDofsLocalWithoutGhosts: " << nDofsLocalWithoutGhosts << ", nDofsGlobal: " << nDofsGlobal << ", nGhostDofs: " << nGhostDofs << ", localIndices: " << localIndices;

  // set values in vector
  // initialize PETSc vector object
  //ierr = VecCreate(this->meshPartition_->mpiCommunicator(), &vectorLocal_[componentNo]); CHKERRQ(ierr);

  // set sparsity type and other options
  ierr = VecSetFromOptions(rhs); CHKERRQ(ierr);
  //ierr = VecGhostGetLocalForm(rhs, &rhsLocal); CHKERRQ(ierr);

  ierr = VecZeroEntries(rhs); CHKERRQ(ierr);

  VecAssemblyBegin(rhs);
  VecAssemblyEnd(rhs);
  if (ownRankNo == 0)
  {
    // set the values
    std::array<double,8> values = {1.0,  1.0,  1.0,  1.0,  -0.333333,  -0.333333,  -0.333333,  -0.333333};
    ierr = VecSetValues(rhs, 8, localIndices.data(), values.data(), INSERT_VALUES);
  }
  else if (ownRankNo == 2)
  {
    // set the values
    std::array<double,8> values = {-0.666667, -0.666667, -0.666667, -0.666667, 2.0, 2.0, 2.0, 2.0};
    ierr = VecSetValues(rhs, 8, localIndices.data(), values.data(), INSERT_VALUES);
  }

  VecAssemblyBegin(rhs);
  VecAssemblyEnd(rhs);

  VecGhostUpdateBegin(rhs,INSERT_VALUES,SCATTER_FORWARD);
  VecGhostUpdateEnd(rhs,INSERT_VALUES,SCATTER_FORWARD);

  // output rhs vector
  std::vector<double> vectorValues(nDofsLocalWithoutGhosts);
  VecGetValues(rhs, nDofsLocalWithoutGhosts, localIndices.data(), vectorValues.data());
  LOG(DEBUG) << "rhs: " << vectorValues;


  Mat mat;

  int nNonZerosDiagonal = 27;
  int nNonZerosOffdiagonal = 27;

  dof_no_t nRowsLocal = nDofsLocalWithoutGhosts;
  dof_no_t nRowsGlobal = nDofsGlobal;

  dof_no_t nColumnsLocal = nDofsLocalWithoutGhosts;
  dof_no_t nColumnsGlobal = nDofsGlobal;

  //ierr = MatCreateAIJ(rankSubset_->mpiCommunicator(), partition.(), partition.(), n, n,
  //                    nNonZerosDiagonal, NULL, nNonZerosOffdiagonal, NULL, &matrix); CHKERRQ(ierr);

  // parallel API
  //ierr = DMSetMatrixPreallocateOnly(this->dm_, PETSC_TRUE); CHKERRQ(ierr);  // do not fill zero entries when DMCreateMatrix is called
  //ierr = DMCreateMatrix(dm_, &matrix_); CHKERRQ(ierr);

  ierr = MatCreate(mpiCommunicator, &mat); CHKERRQ(ierr);
  ierr = MatSetSizes(mat, nRowsLocal, nColumnsLocal,
                      nRowsGlobal, nColumnsGlobal); CHKERRQ(ierr);

  ierr = MatSetType(mat, MATAIJ); CHKERRQ(ierr);

  //ierr = MatSetFromOptions(mat); CHKERRQ(ierr);         // either use MatSetFromOptions or MatSetUp to allocate internal data structures

  // allow additional non-zero entries in the stiffness matrix for UnstructuredDeformable mesh
  //MatSetOption(matrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);

  // dense matrix
  // MATDENSE = "dense" - A matrix type to be used for dense matrices. This matrix type is identical to MATSEQDENSE when constructed with a single process communicator, and MATMPIDENSE otherwise.
  //ierr = MatSetUp(mat); CHKERRQ(ierr);

  // sparse matrix: preallocation of internal data structure
  // http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MATAIJ.html#MATAIJ
  // MATAIJ = "aij" - A matrix type to be used for sparse matrices. This matrix type is identical to MATSEQAIJ when constructed with a single process communicator, and MATMPIAIJ otherwise.
  // As a result, for single process communicators, MatSeqAIJSetPreallocation is supported, and similarly MatMPIAIJSetPreallocation is supported for communicators controlling multiple processes.
  // It is recommended that you call both of the above preallocation routines for simplicity.
  ierr = MatSeqAIJSetPreallocation(mat, nNonZerosDiagonal, NULL); CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(mat, nNonZerosDiagonal, NULL, nNonZerosOffdiagonal, NULL); CHKERRQ(ierr);
  LOG(DEBUG) << "Mat SetPreallocation, nNonZerosDiagonal: " << nNonZerosDiagonal << ", nNonZerosOffdiagonal: " << nNonZerosOffdiagonal;

  // set values
  std::vector<int> rowIndices, columnIndices;
  if (ownRankNo == 0)
  {
    MatSetValue(mat,0,0,1,INSERT_VALUES);
    MatSetValue(mat,1,1,1,INSERT_VALUES);
    MatSetValue(mat,2,2,1,INSERT_VALUES);
    MatSetValue(mat,3,3,1,INSERT_VALUES);

    std::array<int,4> indices = {4,5,6,7};
    std::array<double,16> matValues =
    {
      -0.63, -0.0648, -0.0648, +0.0926,
      -0.0648, -0.63, +0.0926, -0.0648,
      -0.0648, +0.0926, -0.63, -0.0648,
      +0.0926, -0.0648, -0.0648, -0.63
    };
    MatSetValues(mat, 4, indices.data(), 4, indices.data(), matValues.data(), INSERT_VALUES);

    LOG(DEBUG) << " set values at indices " << indices;

    rowIndices = std::vector<int>({0,1,2,3,4,5,6,7});
    columnIndices = std::vector<int>({0,1,2,3,4,5,6,7});


//     _____0_____1_____2_____3_____4_____5_____6_____7______
//   0|    +1
//  +1|          +1
//  +2|                +1
//  +3|                      +1
//  +4|                         -0.63 -0.0648 -0.0648 +0.0926
//  +5|                         -0.0648 -0.63 +0.0926 -0.0648
//  +6|                         -0.0648 +0.0926 -0.63 -0.0648
//  +7|                         +0.0926 -0.0648 -0.0648 -0.63

  }
  else if (ownRankNo == 1)
  {
    rowIndices = std::vector<int>({8,9,10,11});
    columnIndices = std::vector<int>({8,9,10,11});

    std::array<double,16> matValues =
    {
      -0.63, -0.0648, -0.0648, +0.0926,
      -0.0648, -0.63, +0.0926, -0.0648,
      -0.0648, +0.0926, -0.63, -0.0648,
      +0.0926, -0.0648, -0.0648, -0.63
    };
    MatSetValues(mat, 4, rowIndices.data(), 4, columnIndices.data(), matValues.data(), INSERT_VALUES);

    LOG(DEBUG) << " set values at indices " << rowIndices << ", " << columnIndices;
//     _____0_____1_____2_____3______
//   0| -0.63 -0.0648 -0.0648 +0.0926
//  +1| -0.0648 -0.63 +0.0926 -0.0648
//  +2| -0.0648 +0.0926 -0.63 -0.0648
//  +3| +0.0926 -0.0648 -0.0648 -0.63

  }
  else if (ownRankNo == 2)
  {
    rowIndices = std::vector<int>({12,13,14,15,16,17,18,19});
    columnIndices = std::vector<int>({12,13,14,15,16,17,18,19});

    std::array<int,4> matIndices = {12,13,14,15};
    std::array<double,16> matValues =
    {
      -0.63, -0.0648, -0.0648, +0.0926,
      -0.0648, -0.63, +0.0926, -0.0648,
      -0.0648, +0.0926, -0.63, -0.0648,
      +0.0926, -0.0648, -0.0648, -0.63
    };
    MatSetValues(mat, 4, matIndices.data(), 4, matIndices.data(), matValues.data(), INSERT_VALUES);

    MatSetValue(mat,16,16,1,INSERT_VALUES);
    MatSetValue(mat,17,17,1,INSERT_VALUES);
    MatSetValue(mat,18,18,1,INSERT_VALUES);
    MatSetValue(mat,19,19,1,INSERT_VALUES);

    LOG(DEBUG) << " set values at indices " << rowIndices << ", " << columnIndices;
//     _____0_____1_____2_____3_____4_____5_____6_____7_____
//   0|    +1
//  +1|          +1
//  +2|                +1
//  +3|                      +1
//  +4|                         -0.63 -0.0648 -0.0648 +0.0926
//  +5|                         -0.0648 -0.63 +0.0926 -0.0648
//  +6|                         -0.0648 +0.0926 -0.63 -0.0648
//  +7|                         +0.0926 -0.0648 -0.0648 -0.63

  }

  MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY);

  {
    LOG(DEBUG) << "matrix: nRowsLocal: " << nRowsLocal << ", rowIndices: " << rowIndices << ", nColumnsLocal: " << nColumnsLocal << ", columnIndices: " << columnIndices;

    std::vector<double> matrixValues(nDofsLocalWithoutGhosts*nDofsLocalWithoutGhosts);
    // get values in row-major format with global indexing, note, there is no "MatGetValuesLocal"
    MatGetValues(mat, nDofsLocalWithoutGhosts, localIndices.data(), nDofsLocalWithoutGhosts, localIndices.data(), matrixValues.data());

    LOG(DEBUG) << PetscUtility::getStringMatrix(matrixValues, nRowsLocal, nColumnsLocal, nDofsGlobal, nDofsGlobal);
  }

  Vec solution;
  VecDuplicate(rhs, &solution);

  KSP ksp;
  ierr = KSPCreate (mpiCommunicator, &ksp); CHKERRQ(ierr);

  // parse the solver and preconditioner types from settings
  KSPType kspType = KSPGMRES;
  PCType pcType = PCNONE;

  // set solver type
  ierr = KSPSetType(ksp, kspType); CHKERRQ(ierr);

  // set options from command line, this overrides the python config
  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

  // extract preconditioner context
  PC pc;
  ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);

  // set type of preconditioner
  ierr = PCSetType(pc, pcType); CHKERRQ(ierr);

  // set options from command line, this overrides the python config
  ierr = PCSetFromOptions(pc); CHKERRQ(ierr);

  //                                    relative tol,      absolute tol,  diverg tol.,   max_iterations
  ierr = KSPSetTolerances (ksp, 1e-5, PETSC_DEFAULT, PETSC_DEFAULT, 10000); CHKERRQ(ierr);

  ierr = KSPSetOperators(ksp, mat, mat); CHKERRQ(ierr);

  LOG(DEBUG) << "solve...";
  // solve the system
  ierr = KSPSolve(ksp, rhs, solution); CHKERRQ(ierr);

  int numberOfIterations = 0;
  PetscReal residualNorm = 0.0;
  ierr = KSPGetIterationNumber(ksp, &numberOfIterations); CHKERRQ(ierr);
  ierr = KSPGetResidualNorm(ksp, &residualNorm); CHKERRQ(ierr);

  KSPConvergedReason convergedReason;
  ierr = KSPGetConvergedReason(ksp, &convergedReason); CHKERRQ(ierr);

  LOG(INFO) << "Solution obtained in " << numberOfIterations << " iterations, residual norm " << residualNorm
    << ": " << PetscUtility::getStringLinearConvergedReason(convergedReason);


  // output solution vector
  vectorValues.resize(nDofsLocalWithoutGhosts);
  VecGetValues(solution, nDofsLocalWithoutGhosts, localIndices.data(), vectorValues.data());
  LOG(DEBUG) << "solution: " << vectorValues;


  MPI_Barrier(mpiCommunicator);

  LOG(INFO) << "-------------------- END TEST ------------------------";
  return EXIT_SUCCESS;
}

int testPartitioned(int argc, char *argv[])
{

  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  // PetscInitialize(&argc,&argv, NULL, NULL);

  PetscErrorCode ierr;
  //int ownRankNo = 0;
  int ownRankNo = settings.ownRankNo();
  std::vector<int> localIndices;
  //MPI_Comm_rank(mpiCommunicator, &ownRankNo);

  typedef FunctionSpace::FunctionSpace< Mesh::StructuredDeformableOfDimension<3>, BasisFunction::LagrangeOfOrder<1>> FunctionSpaceType;

  /*
   *
  FunctionSpaceDofsNodes(std::shared_ptr<Partition::Manager> partitionManager, const std::vector<Vec3> &nodePositions,
                         const std::array<element_no_t,D> nElementsPerCoordinateDirection, const std::array<int,D> nRanksPerCoordinateDirection);
                         */
  std::shared_ptr<Partition::Manager> partitionManager = settings.partitionManager();
  std::vector<Vec3> nodePositions = {
    Vec3({0,0,0}), Vec3({1,0,0}), Vec3({0,1,0}), Vec3({1,1,0}),
    Vec3({0,0,1}), Vec3({1,0,1}), Vec3({0,1,1}), Vec3({1,1,1}),
    Vec3({0,0,2}), Vec3({1,0,2}), Vec3({0,1,2}), Vec3({1,1,2}),
    Vec3({0,0,3}), Vec3({1,0,3}), Vec3({0,1,3}), Vec3({1,1,3}),
    Vec3({0,0,4}), Vec3({1,0,4}), Vec3({0,1,4}), Vec3({1,1,4})
  };
  std::array<element_no_t,3> nElementsPerCoordinateDirection({1,1,4});
  std::array<int,3> nRanksPerCoordinateDirection({1,1,3});
  std::shared_ptr<FunctionSpaceType> functionSpace
    = std::make_shared<FunctionSpaceType>(partitionManager, nodePositions, nElementsPerCoordinateDirection, nRanksPerCoordinateDirection);

  /*
   * createPartitioningStructuredLocal(std::array<global_no_t,FunctionSpace::dim()> &nElementsGlobal,
                                                                                const std::array<element_no_t,FunctionSpace::dim()> nElementsLocal,
                                                                                const std::array<int,FunctionSpace::dim()> nRanks);
                                                                                */
  std::array<global_no_t,3> nElementsGlobal({1,1,4});
  std::array<int,3> nElementsLocal({1,1,1});

  if (ownRankNo == 0)
    nElementsLocal[2] = 2;
  std::vector<int> rankNos;
  PythonConfig specificSettings;
  std::shared_ptr<Partition::MeshPartition<FunctionSpaceType>> meshPartition = settings.partitionManager()->createPartitioningStructuredLocal<FunctionSpaceType>(
    specificSettings, nElementsGlobal, nElementsLocal, nRanksPerCoordinateDirection, rankNos);

  int nDofsGlobal = meshPartition->nDofsGlobal();
  int nDofsLocalWithoutGhosts = meshPartition->nDofsLocalWithoutGhosts();
  int nColumnsLocal = meshPartition->nDofsLocalWithoutGhosts();
  int nRowsLocal = meshPartition->nDofsLocalWithoutGhosts();

  PartitionedPetscVec<FunctionSpaceType,1> rhs(meshPartition, "rhs");
  //PartitionedPetscVec<FunctionSpaceType,1> rhs(*(functionSpace->geometryField(). partitionedPetscVec()), "rhs");

  // set values of rhs
  if (ownRankNo == 0)
  {
    // set the values
    std::array<double,8> values = {1.0,  1.0,  1.0,  1.0,  -0.333333,  -0.333333,  -0.333333,  -0.333333};
    std::vector<int> localIndices = {0,1,2,3,4,5,6,7,8};
    rhs.setValues(0, 8, localIndices.data(), values.data(), INSERT_VALUES);
    // setValues(int componentNo, PetscInt ni, const PetscInt ix[], const PetscScalar y[], InsertMode iora)
  }
  else if (ownRankNo == 2)
  {
    // set the values
    std::array<double,8> values = {-0.666667, -0.666667, -0.666667, -0.666667, 2.0, 2.0, 2.0, 2.0};
    std::vector<int> localIndices = {0,1,2,3,4,5,6,7,8};
    rhs.setValues(0, 8, localIndices.data(), values.data(), INSERT_VALUES);
  }
  std::stringstream str;
  rhs.output(str);
  LOG(DEBUG) << str.str();

  // create the matrix
  /*
   *
  PartitionedPetscMat(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>> meshPartition,
                      int nComponents, int nNonZerosDiagonal, int nNonZerosOffdiagonal, std::string name);

   */
  PartitionedPetscMat<FunctionSpaceType> mat(meshPartition, 1, 27, 27, "stiffnessMatrix");


  if (ownRankNo == 0)
  {
    localIndices = std::vector<int>({0,1,2,3,4,5,6,7});
    mat.setValue(0,0,1,INSERT_VALUES);
    mat.setValue(1,1,1,INSERT_VALUES);
    mat.setValue(2,2,1,INSERT_VALUES);
    mat.setValue(3,3,1,INSERT_VALUES);

    std::array<int,4> indices = {4,5,6,7};
    std::array<double,16> matValues =
    {
      -0.63, -0.0648, -0.0648, +0.0926,
      -0.0648, -0.63, +0.0926, -0.0648,
      -0.0648, +0.0926, -0.63, -0.0648,
      +0.0926, -0.0648, -0.0648, -0.63
    };
    mat.setValues(4, indices.data(), 4, indices.data(), matValues.data(), INSERT_VALUES);

    LOG(DEBUG) << " set values at indices " << indices;


//     _____0_____1_____2_____3_____4_____5_____6_____7______
//   0|    +1
//  +1|          +1
//  +2|                +1
//  +3|                      +1
//  +4|                         -0.63 -0.0648 -0.0648 +0.0926
//  +5|                         -0.0648 -0.63 +0.0926 -0.0648
//  +6|                         -0.0648 +0.0926 -0.63 -0.0648
//  +7|                         +0.0926 -0.0648 -0.0648 -0.63

  }
  else if (ownRankNo == 1)
  {
    localIndices = std::vector<int>({8,9,10,11});
    std::array<double,16> matValues =
    {
      -0.63, -0.0648, -0.0648, +0.0926,
      -0.0648, -0.63, +0.0926, -0.0648,
      -0.0648, +0.0926, -0.63, -0.0648,
      +0.0926, -0.0648, -0.0648, -0.63
    };
    std::vector<int> localIndices({0,1,2,3});
    mat.setValues(4, localIndices.data(), 4, localIndices.data(), matValues.data(), INSERT_VALUES);

//     _____0_____1_____2_____3______
//   0| -0.63 -0.0648 -0.0648 +0.0926
//  +1| -0.0648 -0.63 +0.0926 -0.0648
//  +2| -0.0648 +0.0926 -0.63 -0.0648
//  +3| +0.0926 -0.0648 -0.0648 -0.63

  }
  else if (ownRankNo == 2)
  {
    localIndices = std::vector<int>({12,13,14,15,16,17,18,19});
    std::array<double,16> matValues =
    {
      -0.63, -0.0648, -0.0648, +0.0926,
      -0.0648, -0.63, +0.0926, -0.0648,
      -0.0648, +0.0926, -0.63, -0.0648,
      +0.0926, -0.0648, -0.0648, -0.63
    };
    std::vector<int> localIndices({0,1,2,3});
    mat.setValues( 4, localIndices.data(), 4, localIndices.data(), matValues.data(), INSERT_VALUES);

    mat.setValue(4,4,1,INSERT_VALUES);
    mat.setValue(5,5,1,INSERT_VALUES);
    mat.setValue(6,6,1,INSERT_VALUES);
    mat.setValue(7,7,1,INSERT_VALUES);

//     _____0_____1_____2_____3_____4_____5_____6_____7_____
//   0|    +1
//  +1|          +1
//  +2|                +1
//  +3|                      +1
//  +4|                         -0.63 -0.0648 -0.0648 +0.0926
//  +5|                         -0.0648 -0.63 +0.0926 -0.0648
//  +6|                         -0.0648 +0.0926 -0.63 -0.0648
//  +7|                         +0.0926 -0.0648 -0.0648 -0.63

  }
  mat.assembly(MAT_FINAL_ASSEMBLY);

  LOG(DEBUG) << mat;


  {
    std::vector<double> matrixValues(nDofsLocalWithoutGhosts*nDofsLocalWithoutGhosts);
    // get values in row-major format with global indexing, note, there is no "MatGetValuesLocal"
    MatGetValues(mat.valuesGlobal(), nDofsLocalWithoutGhosts, localIndices.data(), nDofsLocalWithoutGhosts, localIndices.data(), matrixValues.data());

    LOG(DEBUG) << PetscUtility::getStringMatrix(matrixValues, nRowsLocal, nColumnsLocal, nDofsGlobal, nDofsGlobal);
  }


  Vec solution;
  VecDuplicate(rhs.valuesGlobal(), &solution);
/*
  KSP ksp;
  ierr = KSPCreate (mpiCommunicator, &ksp); CHKERRQ(ierr);

  // parse the solver and preconditioner types from settings
  KSPType kspType = KSPGMRES;
  PCType pcType = PCNONE;

  // set solver type
  ierr = KSPSetType(ksp, kspType); CHKERRQ(ierr);

  // set options from command line, this overrides the python config
  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

  // extract preconditioner context
  PC pc;
  ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);

  // set type of preconditioner
  ierr = PCSetType(pc, pcType); CHKERRQ(ierr);

  // set options from command line, this overrides the python config
  ierr = PCSetFromOptions(pc); CHKERRQ(ierr);

  //                                    relative tol,      absolute tol,  diverg tol.,   max_iterations
  ierr = KSPSetTolerances (ksp, 1e-5, PETSC_DEFAULT, PETSC_DEFAULT, 10000); CHKERRQ(ierr);

  ierr = KSPSetOperators(ksp, mat.valuesGlobal(), mat.valuesGlobal()); CHKERRQ(ierr);

  LOG(DEBUG) << "solve...";
  // solve the system
  ierr = KSPSolve(ksp, rhs.valuesGlobal(), solution); CHKERRQ(ierr);

  int numberOfIterations = 0;
  PetscReal residualNorm = 0.0;
  ierr = KSPGetIterationNumber(ksp, &numberOfIterations); CHKERRQ(ierr);
  ierr = KSPGetResidualNorm(ksp, &residualNorm); CHKERRQ(ierr);

  KSPConvergedReason convergedReason;
  ierr = KSPGetConvergedReason(ksp, &convergedReason); CHKERRQ(ierr);

  LOG(INFO) << "Solution obtained in " << numberOfIterations << " iterations, residual norm " << residualNorm
    << ": " << PetscUtility::getStringLinearConvergedReason(convergedReason);
*/

  // get linear solver context from solver manager
  std::shared_ptr<Solver::Linear> linearSolver = settings.solverManager()->template solver<Solver::Linear>(
    settings.getPythonConfig(), meshPartition->mpiCommunicator());
  std::shared_ptr<KSP> ksp = linearSolver->ksp();
  assert(ksp != nullptr);

  // set matrix used for linear system and preconditioner to ksp context
  ierr = KSPSetOperators(*ksp, mat.valuesGlobal(), mat.valuesGlobal()); CHKERRQ(ierr);

  // non-zero initial values
#if 0
  PetscScalar scalar = 0.5;
  ierr = VecSet(data_.solution()->values(), scalar); CHKERRQ(ierr);
  ierr = KSPSetInitialGuessNonzero(*ksp, PETSC_TRUE); CHKERRQ(ierr);
#endif

  LOG(DEBUG) << "solve...";
  // solve the system
  ierr = KSPSolve(*ksp, rhs.valuesGlobal(), solution); CHKERRQ(ierr);

  int numberOfIterations = 0;
  PetscReal residualNorm = 0.0;
  ierr = KSPGetIterationNumber(*ksp, &numberOfIterations); CHKERRQ(ierr);
  ierr = KSPGetResidualNorm(*ksp, &residualNorm); CHKERRQ(ierr);

  KSPConvergedReason convergedReason;
  ierr = KSPGetConvergedReason(*ksp, &convergedReason); CHKERRQ(ierr);

  LOG(INFO) << "Solution obtained in " << numberOfIterations << " iterations, residual norm " << residualNorm
    << ": " << PetscUtility::getStringLinearConvergedReason(convergedReason);


  // output solution vector
  std::vector<double> vectorValues(nDofsLocalWithoutGhosts);
  VecGetValues(solution, nDofsLocalWithoutGhosts, localIndices.data(), vectorValues.data());
  LOG(DEBUG) << "solution: " << vectorValues;

  LOG(INFO) << "-------------------- END TEST ------------------------";
  return EXIT_SUCCESS;
}

int testFieldVariables(int argc, char *argv[])
{
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  // PetscInitialize(&argc,&argv, NULL, NULL);

  PetscErrorCode ierr;
  //int ownRankNo = 0;
  int ownRankNo = settings.ownRankNo();
  std::vector<int> localIndices;
  //MPI_Comm_rank(mpiCommunicator, &ownRankNo);

  typedef FunctionSpace::FunctionSpace< Mesh::StructuredDeformableOfDimension<3>, BasisFunction::LagrangeOfOrder<1>> FunctionSpaceType;

  /*
   *
  FunctionSpaceDofsNodes(std::shared_ptr<Partition::Manager> partitionManager, const std::vector<Vec3> &nodePositions,
                         const std::array<element_no_t,D> nElementsPerCoordinateDirection, const std::array<int,D> nRanksPerCoordinateDirection);
                         */
  std::shared_ptr<Partition::Manager> partitionManager = settings.partitionManager();


  // create mesh partition
  /*
   * createPartitioningStructuredLocal(std::array<global_no_t,FunctionSpace::dim()> &nElementsGlobal,
                                                                                const std::array<element_no_t,FunctionSpace::dim()> nElementsLocal,
                                                                                const std::array<int,FunctionSpace::dim()> nRanks);
                                                                                */
  std::array<global_no_t,3> nElementsGlobal({1,1,4});
  std::array<int,3> nElementsLocal({1,1,1});
  std::array<int,3> nRanksPerCoordinateDirection({1,1,3});

  if (ownRankNo == 0)
    nElementsLocal[2] = 2;
  std::vector<int> rankNos;
  PythonConfig specificSettings;
  std::shared_ptr<Partition::MeshPartition<FunctionSpaceType>> meshPartition = settings.partitionManager()->createPartitioningStructuredLocal<FunctionSpaceType>(
    specificSettings, nElementsGlobal, nElementsLocal, nRanksPerCoordinateDirection, rankNos);

  // create function space
  std::vector<Vec3> nodePositions = {
    Vec3({0,0,0}), Vec3({1,0,0}), Vec3({0,1,0}), Vec3({1,1,0}),
    Vec3({0,0,1}), Vec3({1,0,1}), Vec3({0,1,1}), Vec3({1,1,1}),
    Vec3({0,0,2}), Vec3({1,0,2}), Vec3({0,1,2}), Vec3({1,1,2}),
    Vec3({0,0,3}), Vec3({1,0,3}), Vec3({0,1,3}), Vec3({1,1,3}),
    Vec3({0,0,4}), Vec3({1,0,4}), Vec3({0,1,4}), Vec3({1,1,4})
  };
  std::array<element_no_t,3> nElementsPerCoordinateDirection({1,1,4});
  std::shared_ptr<FunctionSpaceType> functionSpace
    = settings.meshManager()->createFunctionSpaceWithGivenMeshPartition<FunctionSpaceType>("f", meshPartition, nodePositions, nElementsPerCoordinateDirection, nRanksPerCoordinateDirection);
/*
  std::shared_ptr<FunctionSpaceType> functionSpace
    = std::make_shared<FunctionSpaceType>(partitionManager, nodePositions, nElementsPerCoordinateDirection, nRanksPerCoordinateDirection);
*/

  int nDofsGlobal = meshPartition->nDofsGlobal();
  int nDofsLocalWithoutGhosts = meshPartition->nDofsLocalWithoutGhosts();
  int nColumnsLocal = meshPartition->nDofsLocalWithoutGhosts();
  int nRowsLocal = meshPartition->nDofsLocalWithoutGhosts();

  //std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> rhs
  //= functionSpace->template createFieldVariable<1>("rhs");
  std::vector<std::string> componentNames({"0"});
  /*std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> rhs = std::make_shared<FieldVariable::FieldVariable<FunctionSpaceType,1>>(
    functionSpace, "f", componentNames
  );
  */// works


  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> rhs
    = std::make_shared<FieldVariable::FieldVariable<FunctionSpaceType,1>>(
      functionSpace->geometryField(), "f", componentNames
    );


  rhs->zeroEntries();
  // std::shared_ptr<FunctionSpaceType> functionSpace, std::string name, std::vector<std::string> componentNames, bool isGeometryField=false);
  //std::shared_ptr<FieldVariable<FunctionSpaceType,1>> solution = functionSpace->template createFieldVariable<1>("solution");

  // set values of rhs
  if (ownRankNo == 0)
  {
    // set the values
    std::vector<double> values = {1.0,  1.0,  1.0,  1.0,  -0.333333,  -0.333333,  -0.333333,  -0.333333};
    rhs->setValuesWithoutGhosts(values, INSERT_VALUES);
    // setValues(int componentNo, PetscInt ni, const PetscInt ix[], const PetscScalar y[], InsertMode iora)
  }
  else if (ownRankNo == 2)
  {
    // set the values
    std::vector<double> values = {-0.666667, -0.666667, -0.666667, -0.666667, 2.0, 2.0, 2.0, 2.0};
    rhs->setValuesWithoutGhosts(values, INSERT_VALUES);
  }
  std::stringstream str;
  rhs->output(str);
  LOG(DEBUG) << str.str();

  // create the matrix
  /*
   *
  PartitionedPetscMat(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>> meshPartition,
                      int nComponents, int nNonZerosDiagonal, int nNonZerosOffdiagonal, std::string name);

   */
  PartitionedPetscMat<FunctionSpaceType> mat(meshPartition, 1, 27, 27, "stiffnessMatrix");


  if (ownRankNo == 0)
  {
    localIndices = std::vector<int>({0,1,2,3,4,5,6,7});
    mat.setValue(0,0,1,INSERT_VALUES);
    mat.setValue(1,1,1,INSERT_VALUES);
    mat.setValue(2,2,1,INSERT_VALUES);
    mat.setValue(3,3,1,INSERT_VALUES);

    std::array<int,4> indices = {4,5,6,7};
    std::array<double,16> matValues =
    {
      -0.63, -0.0648, -0.0648, +0.0926,
      -0.0648, -0.63, +0.0926, -0.0648,
      -0.0648, +0.0926, -0.63, -0.0648,
      +0.0926, -0.0648, -0.0648, -0.63
    };
    mat.setValues(4, indices.data(), 4, indices.data(), matValues.data(), INSERT_VALUES);

    LOG(DEBUG) << " set values at indices " << indices;


//     _____0_____1_____2_____3_____4_____5_____6_____7______
//   0|    +1
//  +1|          +1
//  +2|                +1
//  +3|                      +1
//  +4|                         -0.63 -0.0648 -0.0648 +0.0926
//  +5|                         -0.0648 -0.63 +0.0926 -0.0648
//  +6|                         -0.0648 +0.0926 -0.63 -0.0648
//  +7|                         +0.0926 -0.0648 -0.0648 -0.63

  }
  else if (ownRankNo == 1)
  {
    localIndices = std::vector<int>({8,9,10,11});
    std::array<double,16> matValues =
    {
      -0.63, -0.0648, -0.0648, +0.0926,
      -0.0648, -0.63, +0.0926, -0.0648,
      -0.0648, +0.0926, -0.63, -0.0648,
      +0.0926, -0.0648, -0.0648, -0.63
    };
    std::vector<int> localIndices({0,1,2,3});
    mat.setValues(4, localIndices.data(), 4, localIndices.data(), matValues.data(), INSERT_VALUES);

//     _____0_____1_____2_____3______
//   0| -0.63 -0.0648 -0.0648 +0.0926
//  +1| -0.0648 -0.63 +0.0926 -0.0648
//  +2| -0.0648 +0.0926 -0.63 -0.0648
//  +3| +0.0926 -0.0648 -0.0648 -0.63

  }
  else if (ownRankNo == 2)
  {
    localIndices = std::vector<int>({12,13,14,15,16,17,18,19});
    std::array<double,16> matValues =
    {
      -0.63, -0.0648, -0.0648, +0.0926,
      -0.0648, -0.63, +0.0926, -0.0648,
      -0.0648, +0.0926, -0.63, -0.0648,
      +0.0926, -0.0648, -0.0648, -0.63
    };
    std::vector<int> localIndices({0,1,2,3});
    mat.setValues( 4, localIndices.data(), 4, localIndices.data(), matValues.data(), INSERT_VALUES);

    mat.setValue(4,4,1,INSERT_VALUES);
    mat.setValue(5,5,1,INSERT_VALUES);
    mat.setValue(6,6,1,INSERT_VALUES);
    mat.setValue(7,7,1,INSERT_VALUES);

//     _____0_____1_____2_____3_____4_____5_____6_____7_____
//   0|    +1
//  +1|          +1
//  +2|                +1
//  +3|                      +1
//  +4|                         -0.63 -0.0648 -0.0648 +0.0926
//  +5|                         -0.0648 -0.63 +0.0926 -0.0648
//  +6|                         -0.0648 +0.0926 -0.63 -0.0648
//  +7|                         +0.0926 -0.0648 -0.0648 -0.63

  }
  mat.assembly(MAT_FINAL_ASSEMBLY);

  LOG(DEBUG) << mat;


  {
    std::vector<double> matrixValues(nDofsLocalWithoutGhosts*nDofsLocalWithoutGhosts);
    // get values in row-major format with global indexing, note, there is no "MatGetValuesLocal"
    MatGetValues(mat.valuesGlobal(), nDofsLocalWithoutGhosts, localIndices.data(), nDofsLocalWithoutGhosts, localIndices.data(), matrixValues.data());

    LOG(DEBUG) << PetscUtility::getStringMatrix(matrixValues, nRowsLocal, nColumnsLocal, nDofsGlobal, nDofsGlobal);
  }


  Vec solution;
  VecDuplicate(rhs->valuesGlobal(), &solution);
/*
  KSP ksp;
  ierr = KSPCreate (mpiCommunicator, &ksp); CHKERRQ(ierr);

  // parse the solver and preconditioner types from settings
  KSPType kspType = KSPGMRES;
  PCType pcType = PCNONE;

  // set solver type
  ierr = KSPSetType(ksp, kspType); CHKERRQ(ierr);

  // set options from command line, this overrides the python config
  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

  // extract preconditioner context
  PC pc;
  ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);

  // set type of preconditioner
  ierr = PCSetType(pc, pcType); CHKERRQ(ierr);

  // set options from command line, this overrides the python config
  ierr = PCSetFromOptions(pc); CHKERRQ(ierr);

  //                                    relative tol,      absolute tol,  diverg tol.,   max_iterations
  ierr = KSPSetTolerances (ksp, 1e-5, PETSC_DEFAULT, PETSC_DEFAULT, 10000); CHKERRQ(ierr);

  ierr = KSPSetOperators(ksp, mat.valuesGlobal(), mat.valuesGlobal()); CHKERRQ(ierr);

  LOG(DEBUG) << "solve...";
  // solve the system
  ierr = KSPSolve(ksp, rhs->valuesGlobal(), solution); CHKERRQ(ierr);

  int numberOfIterations = 0;
  PetscReal residualNorm = 0.0;
  ierr = KSPGetIterationNumber(ksp, &numberOfIterations); CHKERRQ(ierr);
  ierr = KSPGetResidualNorm(ksp, &residualNorm); CHKERRQ(ierr);

  KSPConvergedReason convergedReason;
  ierr = KSPGetConvergedReason(ksp, &convergedReason); CHKERRQ(ierr);

  LOG(INFO) << "Solution obtained in " << numberOfIterations << " iterations, residual norm " << residualNorm
    << ": " << PetscUtility::getStringLinearConvergedReason(convergedReason);
*/

  // get linear solver context from solver manager
  std::shared_ptr<Solver::Linear> linearSolver = settings.solverManager()->template solver<Solver::Linear>(
    settings.getPythonConfig(), meshPartition->mpiCommunicator());
  std::shared_ptr<KSP> ksp = linearSolver->ksp();
  assert(ksp != nullptr);

  // set matrix used for linear system and preconditioner to ksp context
  ierr = KSPSetOperators(*ksp, mat.valuesGlobal(), mat.valuesGlobal()); CHKERRQ(ierr);

  // non-zero initial values
#if 0
  PetscScalar scalar = 0.5;
  ierr = VecSet(data_.solution()->values(), scalar); CHKERRQ(ierr);
  ierr = KSPSetInitialGuessNonzero(*ksp, PETSC_TRUE); CHKERRQ(ierr);
#endif

  LOG(DEBUG) << "solve...";
  // solve the system
  ierr = KSPSolve(*ksp, rhs->valuesGlobal(), solution); CHKERRQ(ierr);

  int numberOfIterations = 0;
  PetscReal residualNorm = 0.0;
  ierr = KSPGetIterationNumber(*ksp, &numberOfIterations); CHKERRQ(ierr);
  ierr = KSPGetResidualNorm(*ksp, &residualNorm); CHKERRQ(ierr);

  KSPConvergedReason convergedReason;
  ierr = KSPGetConvergedReason(*ksp, &convergedReason); CHKERRQ(ierr);

  LOG(INFO) << "Solution obtained in " << numberOfIterations << " iterations, residual norm " << residualNorm
    << ": " << PetscUtility::getStringLinearConvergedReason(convergedReason);


  // output solution vector
  std::vector<double> vectorValues(nDofsLocalWithoutGhosts);
  VecGetValues(solution, nDofsLocalWithoutGhosts, localIndices.data(), vectorValues.data());
  LOG(DEBUG) << "solution: " << vectorValues;

  LOG(INFO) << "-------------------- END TEST ------------------------";
  return EXIT_SUCCESS;
}

int main(int argc, char *argv[])
{
  return testFieldVariables(argc, argv);
  //return testPartitioned(argc, argv);
}
