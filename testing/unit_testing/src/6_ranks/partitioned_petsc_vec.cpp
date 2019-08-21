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

TEST(PartitionedPetscVecTest, Test)
{
  // explicit functionSpace with node positions
  std::string pythonConfig = R"(
import numpy as np
# Laplace 1D
config = {
  "Meshes" : {
    "testMesh": {
      "nElements": [6,3],
      "physicalExtent": [1.0,1.0],
    }
  },
  "FiniteElementMethod" : {
    "relativeTolerance": 1e-15,
    "meshName": "testMesh",
  },
  "dirichletBoundaryConditions": {6:[1.0,np.nan], 2:[np.nan,2.0], 28:[1.0,np.nan], -4:[2.0,1.0], -12:[3,4], 7: [np.nan,5]},
}
)";

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
  else if (settings.ownRankNo() == 1)
  {
    vec0->setValue(0, 3, 105);
    vec0->setValue(1, 3, 115);
    vec0->setValue(0, 0, 123);
    vec0->setValue(1, 0, 113);
  }
  else if (settings.ownRankNo() == 2)
  {
    vec0->setValue(0, 2, 114);
    vec0->setValue(1, 2, 104);
    vec0->setValue(0, 1, 133);
    vec0->setValue(1, 0, 114);
  }
  else if (settings.ownRankNo() == 3)
  {
    vec0->setValue(0, 3, 142);
    vec0->setValue(1, 3, 109);
    vec0->setValue(0, 0, 132);
  }
  else if (settings.ownRankNo() == 4)
  {
    vec0->setValue(0, 1, 123);
    vec0->setValue(1, 1, 128);
    vec0->setValue(0, 3, 113);
    vec0->setValue(1, 5, 134);
    vec0->setValue(0, 5, 5);
  }
  else if (settings.ownRankNo() == 5)
  {
    vec0->setValue(0, 2, 103);
    vec0->setValue(1, 4, 158);
    vec0->setValue(0, 4, 113);
    vec0->setValue(1, 2, 162);
  }

  vec0->finishGhostManipulation();
  vec0->startGhostManipulation();


  std::vector<int> indices({0, 1, 2, 3});
  std::vector<double> values(4);
  vec0->getValues(0, 4, indices.data(), values.data());

  LOG(DEBUG) << "values component 0: " << values;
  if (settings.ownRankNo() == 0)
  {
    ASSERT_EQ(values[0], 100);
    ASSERT_EQ(values[1], 101);
    ASSERT_EQ(values[2], 102);
    ASSERT_EQ(values[3], 4);
  }
  else if (settings.ownRankNo() == 1)
  {
    ASSERT_EQ(values[0], 123);
    ASSERT_EQ(values[1], 0);
    ASSERT_EQ(values[2], 0);
    ASSERT_EQ(values[3], 105);
  }
  else if (settings.ownRankNo() == 2)
  {
    ASSERT_EQ(values[0], 0);
    ASSERT_EQ(values[1], 133);
    ASSERT_EQ(values[2], 1);
    ASSERT_EQ(values[3], 0);
  }
  else if (settings.ownRankNo() == 3)
  {
    ASSERT_EQ(values[0], 132);
    ASSERT_EQ(values[1], 0);
    ASSERT_EQ(values[2], 0);
    ASSERT_EQ(values[3], 142);
  }

  else if (settings.ownRankNo() == 4)
  {
    ASSERT_EQ(values[0], 3);
    ASSERT_EQ(values[1], 123);
    ASSERT_EQ(values[2], 0);
    ASSERT_EQ(values[3], 2);
  }

  else if (settings.ownRankNo() == 5)
  {
    ASSERT_EQ(values[0], 0);
    ASSERT_EQ(values[1], 0);
    ASSERT_EQ(values[2], 103);
    ASSERT_EQ(values[3], 5);
  }


  vec0->getValues(1, 4, indices.data(), values.data());

  LOG(DEBUG) << "values component 1: " << values;
  if (settings.ownRankNo() == 0)
  {
    ASSERT_EQ(values[0], 110);
    ASSERT_EQ(values[1], 111);
    ASSERT_EQ(values[2], 5);
    ASSERT_EQ(values[3], 0);
  }
  else if (settings.ownRankNo() == 1)
  {
    ASSERT_EQ(values[0], 2);
    ASSERT_EQ(values[1], 0);
    ASSERT_EQ(values[2], 0);
    ASSERT_EQ(values[3], 115);
  }
  else if (settings.ownRankNo() == 2)
  {
    ASSERT_EQ(values[0], 114);
    ASSERT_EQ(values[1], 0);
    ASSERT_EQ(values[2], 104);
    ASSERT_EQ(values[3], 0);
  }

  else if (settings.ownRankNo() == 3)
  {
    ASSERT_EQ(values[0], 0);
    ASSERT_EQ(values[1], 0);
    ASSERT_EQ(values[2], 0);
    ASSERT_EQ(values[3], 109);
  }

  else if (settings.ownRankNo() == 4)
  {
    ASSERT_EQ(values[0], 4);
    ASSERT_EQ(values[1], 128);
    ASSERT_EQ(values[2], 0);
    ASSERT_EQ(values[3], 1);
  }
  else if (settings.ownRankNo() == 5)
  {
    ASSERT_EQ(values[0], 0);
    ASSERT_EQ(values[1], 0);
    ASSERT_EQ(values[2], 162);
    ASSERT_EQ(values[3], 134);
  }


  nFails += ::testing::Test::HasFailure();
}
