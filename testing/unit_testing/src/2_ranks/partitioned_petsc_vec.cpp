#include <Python.h>  // this has to be the first included header

#include <iostream>
#include <cstdlib>
#include <fstream>

#include "gtest/gtest.h"
#include "arg.h"
#include "opendihu.h"
#include "../utility.h"
#include "partition/partitioned_petsc_vec/partitioned_petsc_vec_with_dirichlet_bc.h"
#include "spatial_discretization/boundary_conditions/dirichlet_boundary_conditions.h"

TEST(PartitionedPetscVecTest, TwoElementsBCGlobal)
{
  // explicit functionSpace with node positions
  std::string pythonConfig = R"(
import numpy as np
# Laplace 1D
config = {
  "Meshes" : {
    "testMesh": {
      "nElements": [2,1],
      "physicalExtent": [1.0,1.0],
    }
  },
  "FiniteElementMethod" : {
    "relativeTolerance": 1e-15,
    "meshName": "testMesh",
  },
  "dirichletBoundaryConditions": {0:[1.0,np.nan], 1:[np.nan,2.0], 5:[1.0,np.nan]},
}
)";

  // r0 r1
  //  _ _
  // |_|_|

  // dirichlet values
  // c0
  // +-+-1
  // 1-+-+
  //
  // c1
  // +-+-+
  // +-2-+

  // non-BC indexing
  // c0
  // 0-5-+
  // +-3-4
  //
  // c1
  // 2-7-8
  // 1-+-6


  DihuContext settings(argc, argv, pythonConfig);

  SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<2>,
    BasisFunction::LagrangeOfOrder<1>,
    Quadrature::Gauss<2>,
    Equation::Static::Laplace
  > finiteElementMethod(settings);

  typedef Mesh::StructuredDeformableOfDimension<2> MeshType;
  typedef BasisFunction::LagrangeOfOrder<1> BasisFunctionType;

  typedef FunctionSpace::FunctionSpace<MeshType,BasisFunctionType> FunctionSpaceType;

  std::shared_ptr<FunctionSpaceType> functionSpace = finiteElementMethod.functionSpace();
  functionSpace->initialize();

  const int nComponents = 2;
  std::shared_ptr<SpatialDiscretization::DirichletBoundaryConditions<FunctionSpaceType, nComponents>> dirichletBoundaryConditions
    = std::make_shared<SpatialDiscretization::DirichletBoundaryConditions<FunctionSpaceType, nComponents>>(settings);
  dirichletBoundaryConditions->initialize(settings.getPythonConfig(), functionSpace, "dirichletBoundaryConditions");

  std::shared_ptr<PartitionedPetscVecWithDirichletBc<FunctionSpaceType,2>> vec0 = std::make_shared<PartitionedPetscVecWithDirichletBc<FunctionSpaceType,2>>(
    functionSpace->meshPartition(), dirichletBoundaryConditions, "test");

  vec0->zeroEntries();

  if (settings.ownRankNo() == 0)
  {
    vec0->setValue(0, 0, 100);
    vec0->setValue(1, 0, 110);
    vec0->setValue(0, 1, 101);
    vec0->setValue(1, 1, 111);
    vec0->setValue(0, 2, 102);
    vec0->setValue(1, 2, 112);
    vec0->setValue(0, 3, 4);
  }
  else
  {
    vec0->setValue(0, 3, 105);
    vec0->setValue(1, 3, 115);
    vec0->setValue(0, 0, 100);
    vec0->setValue(1, 0, 110);
    vec0->setValue(0, 2, 5);
  }

  vec0->finishGhostManipulation();
  vec0->startGhostManipulation();

  // set values
  // c0
  // 101-9-1
  // 1-202-+
  //
  // c1
  // 111-+-115
  // 110-2-+


  std::vector<int> indices({0, 1, 2, 3});
  std::vector<double> values(4);
  vec0->getValues(0, 4, indices.data(), values.data());

  LOG(DEBUG) << "values component 0: " << values;
  if (settings.ownRankNo() == 0)
  {
    ASSERT_EQ(values[0], 1);
    ASSERT_EQ(values[1], 101);
    ASSERT_EQ(values[2], 202);
    ASSERT_EQ(values[3], 9);
  }
  else
  {
    ASSERT_EQ(values[0], 202);
    ASSERT_EQ(values[1], 0);
    ASSERT_EQ(values[2], 9);
    ASSERT_EQ(values[3], 1);
  }

  vec0->getValues(1, 4, indices.data(), values.data());

  LOG(DEBUG) << "values component 1: " << values;
  if (settings.ownRankNo() == 0)
  {
    ASSERT_EQ(values[0], 110);
    ASSERT_EQ(values[1], 111);
    ASSERT_EQ(values[2], 2);
    ASSERT_EQ(values[3], 0);
  }
  else
  {
    ASSERT_EQ(values[0], 2);
    ASSERT_EQ(values[1], 0);
    ASSERT_EQ(values[2], 0);
    ASSERT_EQ(values[3], 115);
  }

  nFails += ::testing::Test::HasFailure();
}

TEST(PartitionedPetscVecTest, TwoElementsBCLocal)
{
  // explicit functionSpace with node positions
  std::string pythonConfig = R"(

import sys
rank_no = (int)(sys.argv[-2])
print("rank_no: {}".format(rank_no))

import numpy as np
# Laplace 1D
config = {
  "Meshes" : {
    "testMesh": {
      "nElements": [2,1],
      "physicalExtent": [1.0,1.0],
    }
  },
  "FiniteElementMethod" : {
    "relativeTolerance": 1e-15,
    "meshName": "testMesh",
  },
  "dirichletBoundaryConditions": {0:[1.0,np.nan]} if (rank_no == 0) else {0:[np.nan,2.0], 3:[1.0,np.nan]},
  "inputMeshIsGlobal": False
}
)";

  // r0 r1
  //  _ _
  // |_|_|

  // dirichlet values
  // c0
  // +-+-1
  // 1-+-+
  //
  // c1
  // +-+-+
  // +-2-+

  // non-BC indexing
  // c0
  // 0-5-+
  // +-3-4
  //
  // c1
  // 2-7-8
  // 1-+-6


  DihuContext settings(argc, argv, pythonConfig);

  SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<2>,
    BasisFunction::LagrangeOfOrder<1>,
    Quadrature::Gauss<2>,
    Equation::Static::Laplace
  > finiteElementMethod(settings);

  typedef Mesh::StructuredDeformableOfDimension<2> MeshType;
  typedef BasisFunction::LagrangeOfOrder<1> BasisFunctionType;

  typedef FunctionSpace::FunctionSpace<MeshType,BasisFunctionType> FunctionSpaceType;

  std::shared_ptr<FunctionSpaceType> functionSpace = finiteElementMethod.functionSpace();
  functionSpace->initialize();

  const int nComponents = 2;
  std::shared_ptr<SpatialDiscretization::DirichletBoundaryConditions<FunctionSpaceType, nComponents>> dirichletBoundaryConditions
    = std::make_shared<SpatialDiscretization::DirichletBoundaryConditions<FunctionSpaceType, nComponents>>(settings);
  dirichletBoundaryConditions->initialize(settings.getPythonConfig(), functionSpace, "dirichletBoundaryConditions");

  std::shared_ptr<PartitionedPetscVecWithDirichletBc<FunctionSpaceType,2>> vec0 = std::make_shared<PartitionedPetscVecWithDirichletBc<FunctionSpaceType,2>>(
    functionSpace->meshPartition(), dirichletBoundaryConditions, "test");

  vec0->zeroEntries();

  if (settings.ownRankNo() == 0)
  {
    std::vector<int> indices({0,1,3});
    std::vector<double> values({100, 101, 4});
    vec0->setValue(1, 0, 110);
    vec0->setValues(0, indices.size(), indices.data(), values.data());
    vec0->setValue(1, 1, 111);
    vec0->setValue(0, 2, 102);
    vec0->setValue(1, 2, 112);
  }
  else
  {
    vec0->setValue(0, 3, 105);
    vec0->setValue(0, 0, 100);
    vec0->setValue(1, 0, 110);
    vec0->setValue(0, 2, 5);
  }

  vec0->finishGhostManipulation();
  vec0->startGhostManipulation();
  vec0->zeroGhostBuffer();

  if (settings.ownRankNo() == 1)
  {
    vec0->setValue(1, 3, 15, ADD_VALUES);
    vec0->setValue(1, 3, 100, ADD_VALUES);
  }
  vec0->finishGhostManipulation();
  vec0->startGhostManipulation();

  // set values
  // c0
  // 101-9-1
  // 1-202-+
  //
  // c1
  // 111-+-115
  // 110-2-+


  std::vector<int> indices({0, 1, 2, 3});
  std::vector<double> values(4);
  vec0->getValues(0, 4, indices.data(), values.data());

  LOG(DEBUG) << "values component 0: " << values;
  if (settings.ownRankNo() == 0)
  {
    ASSERT_EQ(values[0], 1);
    ASSERT_EQ(values[1], 101);
    ASSERT_EQ(values[2], 202);
    ASSERT_EQ(values[3], 9);
  }
  else
  {
    ASSERT_EQ(values[0], 202);
    ASSERT_EQ(values[1], 0);
    ASSERT_EQ(values[2], 9);
    ASSERT_EQ(values[3], 1);
  }

  vec0->getValues(1, 4, indices.data(), values.data());

  LOG(DEBUG) << "values component 1: " << values;
  if (settings.ownRankNo() == 0)
  {
    ASSERT_EQ(values[0], 110);
    ASSERT_EQ(values[1], 111);
    ASSERT_EQ(values[2], 2);
    ASSERT_EQ(values[3], 0);
  }
  else
  {
    ASSERT_EQ(values[0], 2);
    ASSERT_EQ(values[1], 0);
    ASSERT_EQ(values[2], 0);
    ASSERT_EQ(values[3], 115);
  }

  nFails += ::testing::Test::HasFailure();
}
