#include <Python.h>  // this has to be the first included header

#include <iostream>
#include <cstdlib>
#include <fstream>

#include "gtest/gtest.h"
#include "arg.h"
#include "opendihu.h"
#include "../utility.h"


TEST(NumberingsTest, DofNodeNumberingDeformableLinear1D)
{
  // methods that are tested in this test case:
  // global_no_t meshPartition->getElementNoGlobalNatural(element_no_t elementNoLocal)
  // std::array<int,D> meshPartition_->getNodeNoGlobalCoordinates(node_no_t nodeNoLocal)
  // global_no_t functionSpace->getNodeNoGlobalNatural(global_no_t elementNoGlobalNatural, int nodeIndex) const
  // node_no_t functionSpace->getNeighbourNodeNoLocal(node_no_t nodeNoLocal, Mesh::face_t direction)
  // node_no_t functionSpace->getNodeNo(element_no_t elementNo, int nodeIndex)
  // meshPartition->getNodeNoGlobalNatural(coordinatesGlobal)

  // run serial problem
  std::string pythonConfig = R"(

nx = 7
config = {
      "FiniteElementMethod": {
        "nElements": nx,
      }
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  typedef SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<1>,
    BasisFunction::LagrangeOfOrder<1>,
    Quadrature::Gauss<1>,
    Equation::None
  > ProblemType;
  ProblemType problem(settings);

  problem.initialize();

  LOG(INFO) << " ===================== run dof no tests ===================== ";

  typedef typename ProblemType::FunctionSpace FunctionSpace;
  const int D = FunctionSpace::dim();

  std::shared_ptr<FunctionSpace> functionSpace = problem.data().functionSpace();
  std::shared_ptr<Partition::MeshPartition<FunctionSpace>> meshPartition = functionSpace->meshPartition();

  // loop over local elements with ghosts
  for (element_no_t elementNoLocal = 0; elementNoLocal < functionSpace->nElementsLocal(); elementNoLocal++)
  {
    global_no_t elementNoGlobalNatural = meshPartition->getElementNoGlobalNatural(elementNoLocal);

    LOG(INFO) << "";
    LOG(INFO) << "elementNoLocal = " << elementNoLocal << ", elementNoGlobalNatural = " << elementNoGlobalNatural;

    // loop over local nodes, with ghosts
    for (int nodeIndex = 0; nodeIndex < functionSpace->nNodesPerElement(); nodeIndex++)
    {
      node_no_t nodeNoLocal = functionSpace->getNodeNo(elementNoLocal, nodeIndex);

      LOG(INFO) << "elementNoLocal = " << elementNoLocal << ", nodeIndex = " << nodeIndex << ", nodeNoLocal = " << nodeNoLocal;

      std::array<int,D> coordinatesGlobal = meshPartition->getNodeNoGlobalCoordinates(nodeNoLocal);

      LOG(INFO) << "coordinatesGlobal: " << coordinatesGlobal;

      global_no_t nodeNoGlobalNaturalFromCoordinates = meshPartition->getNodeNoGlobalNatural(coordinatesGlobal);
      global_no_t nodeNoGlobalNatural = functionSpace->getNodeNoGlobalNatural(elementNoGlobalNatural, nodeIndex);

      LOG(INFO) << "nodeNoGlobalNaturalFromCoordinates: " << nodeNoGlobalNaturalFromCoordinates << ", nodeNoGlobalNatural: " << nodeNoGlobalNatural;

      ASSERT_EQ(nodeNoGlobalNaturalFromCoordinates, nodeNoGlobalNatural);
    }
  }

  LOG(INFO) << " ===================== run neighbour ===================== ";
  // validate getNeighbourNodeNoLocal
  std::vector<std::pair<std::pair<Mesh::face_t,Mesh::face_t>,std::array<int,D>>> complementFaces = {
    {{Mesh::face_t::face0Plus,  Mesh::face_t::face0Minus}, {1}},
    {{Mesh::face_t::face0Minus, Mesh::face_t::face0Plus}, {-1}}
  };

  for (int i = 0; i < complementFaces.size(); i++)
  {
    Mesh::face_t face0 = complementFaces[i].first.first;
    Mesh::face_t face1 = complementFaces[i].first.second;
    std::array<int,D> expectedDifference = complementFaces[i].second;

    // loop over local nodes without ghosts
    for (node_no_t nodeNoLocal = 0; nodeNoLocal < functionSpace->nNodesLocalWithoutGhosts(); nodeNoLocal++)
    {
      node_no_t neighbourNodeNo = functionSpace->getNeighbourNodeNoLocal(nodeNoLocal, face0);
      LOG(INFO) << "nodeNoLocal = " << nodeNoLocal << ", face " << face0 << ", neighbour: " << neighbourNodeNo;

      if (neighbourNodeNo != -1)
      {
        // compare difference between global coordinates of nodes
        std::array<int,D> coordinatesGlobalNeighbour = meshPartition->getNodeNoGlobalCoordinates(neighbourNodeNo);
        std::array<int,D> coordinatesGlobalOrigin = meshPartition->getNodeNoGlobalCoordinates(nodeNoLocal);


        std::array<int,D> realDifference;
        for(int j = 0; j < D; j++)
        {
          realDifference[j] = coordinatesGlobalNeighbour[j] - coordinatesGlobalOrigin[j];
        }

        LOG(INFO) << "coordinates origin: " << coordinatesGlobalOrigin << ", neighbour: " << coordinatesGlobalNeighbour
          << ", difference: " << realDifference << ", expected: " << expectedDifference;

        for(int j = 0; j < D; j++)
        {
          EXPECT_EQ(expectedDifference[j], realDifference[j]);
        }
      }

      // if the neighbour exists and is a non-ghost dof
      if (neighbourNodeNo != -1 && neighbourNodeNo < functionSpace->nNodesLocalWithoutGhosts())
      {
        node_no_t originNode = functionSpace->getNeighbourNodeNoLocal(neighbourNodeNo, face1);

        LOG(INFO) << "  originNode: " << originNode << " face " << face1;
        ASSERT_EQ(nodeNoLocal, originNode);
      }
    }
  }
}


TEST(NumberingsTest, DofNodeNumberingDeformableQuadratic1D)
{
  // methods that are tested in this test case:
  // global_no_t meshPartition->getElementNoGlobalNatural(element_no_t elementNoLocal)
  // std::array<int,D> meshPartition_->getNodeNoGlobalCoordinates(node_no_t nodeNoLocal)
  // global_no_t functionSpace->getNodeNoGlobalNatural(global_no_t elementNoGlobalNatural, int nodeIndex) const
  // node_no_t functionSpace->getNeighbourNodeNoLocal(node_no_t nodeNoLocal, Mesh::face_t direction)
  // node_no_t functionSpace->getNodeNo(element_no_t elementNo, int nodeIndex)
  // meshPartition->getNodeNoGlobalNatural(coordinatesGlobal)

  // run serial problem
  std::string pythonConfig = R"(

nx = 7
config = {
      "FiniteElementMethod": {
        "nElements": nx,
      }
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  typedef SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<1>,
    BasisFunction::LagrangeOfOrder<2>,
    Quadrature::Gauss<1>,
    Equation::None
  > ProblemType;
  ProblemType problem(settings);

  problem.initialize();

  LOG(INFO) << " ===================== run dof no tests ===================== ";

  typedef typename ProblemType::FunctionSpace FunctionSpace;
  const int D = FunctionSpace::dim();

  std::shared_ptr<FunctionSpace> functionSpace = problem.data().functionSpace();
  std::shared_ptr<Partition::MeshPartition<FunctionSpace>> meshPartition = functionSpace->meshPartition();

  // loop over local elements with ghosts
  for (element_no_t elementNoLocal = 0; elementNoLocal < functionSpace->nElementsLocal(); elementNoLocal++)
  {
    global_no_t elementNoGlobalNatural = meshPartition->getElementNoGlobalNatural(elementNoLocal);

    LOG(INFO) << "";
    LOG(INFO) << "elementNoLocal = " << elementNoLocal << ", elementNoGlobalNatural = " << elementNoGlobalNatural;

    // loop over local nodes, with ghosts
    for (int nodeIndex = 0; nodeIndex < functionSpace->nNodesPerElement(); nodeIndex++)
    {
      node_no_t nodeNoLocal = functionSpace->getNodeNo(elementNoLocal, nodeIndex);

      LOG(INFO) << "elementNoLocal = " << elementNoLocal << ", nodeIndex = " << nodeIndex << ", nodeNoLocal = " << nodeNoLocal;

      std::array<int,D> coordinatesGlobal = meshPartition->getNodeNoGlobalCoordinates(nodeNoLocal);

      LOG(INFO) << "coordinatesGlobal: " << coordinatesGlobal;

      global_no_t nodeNoGlobalNaturalFromCoordinates = meshPartition->getNodeNoGlobalNatural(coordinatesGlobal);
      global_no_t nodeNoGlobalNatural = functionSpace->getNodeNoGlobalNatural(elementNoGlobalNatural, nodeIndex);

      LOG(INFO) << "nodeNoGlobalNaturalFromCoordinates: " << nodeNoGlobalNaturalFromCoordinates << ", nodeNoGlobalNatural: " << nodeNoGlobalNatural;

      ASSERT_EQ(nodeNoGlobalNaturalFromCoordinates, nodeNoGlobalNatural);
    }
  }

  LOG(INFO) << " ===================== run neighbour ===================== ";
  // validate getNeighbourNodeNoLocal
  std::vector<std::pair<std::pair<Mesh::face_t,Mesh::face_t>,std::array<int,D>>> complementFaces = {
    {{Mesh::face_t::face0Plus,  Mesh::face_t::face0Minus}, {1}},
    {{Mesh::face_t::face0Minus, Mesh::face_t::face0Plus}, {-1}},
  };

  for (int i = 0; i < complementFaces.size(); i++)
  {
    Mesh::face_t face0 = complementFaces[i].first.first;
    Mesh::face_t face1 = complementFaces[i].first.second;
    std::array<int,D> expectedDifference = complementFaces[i].second;

    // loop over local nodes without ghosts
    for (node_no_t nodeNoLocal = 0; nodeNoLocal < functionSpace->nNodesLocalWithoutGhosts(); nodeNoLocal++)
    {
      node_no_t neighbourNodeNo = functionSpace->getNeighbourNodeNoLocal(nodeNoLocal, face0);
      LOG(INFO) << "nodeNoLocal = " << nodeNoLocal << ", face " << face0 << ", neighbour: " << neighbourNodeNo;

      if (neighbourNodeNo != -1)
      {
        // compare difference between global coordinates of nodes
        std::array<int,D> coordinatesGlobalNeighbour = meshPartition->getNodeNoGlobalCoordinates(neighbourNodeNo);
        std::array<int,D> coordinatesGlobalOrigin = meshPartition->getNodeNoGlobalCoordinates(nodeNoLocal);


        std::array<int,D> realDifference;
        for(int j = 0; j < D; j++)
        {
          realDifference[j] = coordinatesGlobalNeighbour[j] - coordinatesGlobalOrigin[j];
        }

        LOG(INFO) << "coordinates origin: " << coordinatesGlobalOrigin << ", neighbour: " << coordinatesGlobalNeighbour
          << ", difference: " << realDifference << ", expected: " << expectedDifference;

        for(int j = 0; j < D; j++)
        {
          EXPECT_EQ(expectedDifference[j], realDifference[j]);
        }
      }

      // if the neighbour exists and is a non-ghost dof
      if (neighbourNodeNo != -1 && neighbourNodeNo < functionSpace->nNodesLocalWithoutGhosts())
      {
        node_no_t originNode = functionSpace->getNeighbourNodeNoLocal(neighbourNodeNo, face1);

        LOG(INFO) << "  originNode: " << originNode << " face " << face1;
        ASSERT_EQ(nodeNoLocal, originNode);
      }
    }
  }
}


TEST(NumberingsTest, DofNodeNumberingDeformableHermite1D)
{
  // methods that are tested in this test case:
  // global_no_t meshPartition->getElementNoGlobalNatural(element_no_t elementNoLocal)
  // std::array<int,D> meshPartition_->getNodeNoGlobalCoordinates(node_no_t nodeNoLocal)
  // global_no_t functionSpace->getNodeNoGlobalNatural(global_no_t elementNoGlobalNatural, int nodeIndex) const
  // node_no_t functionSpace->getNeighbourNodeNoLocal(node_no_t nodeNoLocal, Mesh::face_t direction)
  // node_no_t functionSpace->getNodeNo(element_no_t elementNo, int nodeIndex)
  // meshPartition->getNodeNoGlobalNatural(coordinatesGlobal)

  // run serial problem
  std::string pythonConfig = R"(

nx = 7
config = {
      "FiniteElementMethod": {
        "nElements": nx,
      }
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  typedef SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<1>,
    BasisFunction::Hermite,
    Quadrature::Gauss<1>,
    Equation::None
  > ProblemType;
  ProblemType problem(settings);

  problem.initialize();

  LOG(INFO) << " ===================== run dof no tests ===================== ";

  typedef typename ProblemType::FunctionSpace FunctionSpace;
  const int D = FunctionSpace::dim();

  std::shared_ptr<FunctionSpace> functionSpace = problem.data().functionSpace();
  std::shared_ptr<Partition::MeshPartition<FunctionSpace>> meshPartition = functionSpace->meshPartition();

  // loop over local elements with ghosts
  for (element_no_t elementNoLocal = 0; elementNoLocal < functionSpace->nElementsLocal(); elementNoLocal++)
  {
    global_no_t elementNoGlobalNatural = meshPartition->getElementNoGlobalNatural(elementNoLocal);

    LOG(INFO) << "";
    LOG(INFO) << "elementNoLocal = " << elementNoLocal << ", elementNoGlobalNatural = " << elementNoGlobalNatural;

    // loop over local nodes, with ghosts
    for (int nodeIndex = 0; nodeIndex < functionSpace->nNodesPerElement(); nodeIndex++)
    {
      node_no_t nodeNoLocal = functionSpace->getNodeNo(elementNoLocal, nodeIndex);

      LOG(INFO) << "elementNoLocal = " << elementNoLocal << ", nodeIndex = " << nodeIndex << ", nodeNoLocal = " << nodeNoLocal;

      std::array<int,D> coordinatesGlobal = meshPartition->getNodeNoGlobalCoordinates(nodeNoLocal);

      LOG(INFO) << "coordinatesGlobal: " << coordinatesGlobal;

      global_no_t nodeNoGlobalNaturalFromCoordinates = meshPartition->getNodeNoGlobalNatural(coordinatesGlobal);
      global_no_t nodeNoGlobalNatural = functionSpace->getNodeNoGlobalNatural(elementNoGlobalNatural, nodeIndex);

      LOG(INFO) << "nodeNoGlobalNaturalFromCoordinates: " << nodeNoGlobalNaturalFromCoordinates << ", nodeNoGlobalNatural: " << nodeNoGlobalNatural;

      ASSERT_EQ(nodeNoGlobalNaturalFromCoordinates, nodeNoGlobalNatural);
    }
  }
  LOG(INFO) << " ===================== run neighbour ===================== ";
  // validate getNeighbourNodeNoLocal
  std::vector<std::pair<std::pair<Mesh::face_t,Mesh::face_t>,std::array<int,D>>> complementFaces = {
    {{Mesh::face_t::face0Plus,  Mesh::face_t::face0Minus}, {1}},
    {{Mesh::face_t::face0Minus, Mesh::face_t::face0Plus}, {-1}},
  };

  for (int i = 0; i < complementFaces.size(); i++)
  {
    Mesh::face_t face0 = complementFaces[i].first.first;
    Mesh::face_t face1 = complementFaces[i].first.second;
    std::array<int,D> expectedDifference = complementFaces[i].second;

    // loop over local nodes without ghosts
    for (node_no_t nodeNoLocal = 0; nodeNoLocal < functionSpace->nNodesLocalWithoutGhosts(); nodeNoLocal++)
    {
      node_no_t neighbourNodeNo = functionSpace->getNeighbourNodeNoLocal(nodeNoLocal, face0);
      LOG(INFO) << "nodeNoLocal = " << nodeNoLocal << ", face " << face0 << ", neighbour: " << neighbourNodeNo;

      if (neighbourNodeNo != -1)
      {
        // compare difference between global coordinates of nodes
        std::array<int,D> coordinatesGlobalNeighbour = meshPartition->getNodeNoGlobalCoordinates(neighbourNodeNo);
        std::array<int,D> coordinatesGlobalOrigin = meshPartition->getNodeNoGlobalCoordinates(nodeNoLocal);


        std::array<int,D> realDifference;
        for(int j = 0; j < D; j++)
        {
          realDifference[j] = coordinatesGlobalNeighbour[j] - coordinatesGlobalOrigin[j];
        }

        LOG(INFO) << "coordinates origin: " << coordinatesGlobalOrigin << ", neighbour: " << coordinatesGlobalNeighbour
          << ", difference: " << realDifference << ", expected: " << expectedDifference;

        for(int j = 0; j < D; j++)
        {
          EXPECT_EQ(expectedDifference[j], realDifference[j]);
        }
      }

      // if the neighbour exists and is a non-ghost dof
      if (neighbourNodeNo != -1 && neighbourNodeNo < functionSpace->nNodesLocalWithoutGhosts())
      {
        node_no_t originNode = functionSpace->getNeighbourNodeNoLocal(neighbourNodeNo, face1);

        LOG(INFO) << "  originNode: " << originNode << " face " << face1;
        ASSERT_EQ(nodeNoLocal, originNode);
      }
    }
  }
}


// ---------- 2D ---------------

TEST(NumberingsTest, DofNodeNumberingDeformableLinear2D)
{
  // methods that are tested in this test case:
  // global_no_t meshPartition->getElementNoGlobalNatural(element_no_t elementNoLocal)
  // std::array<int,D> meshPartition_->getNodeNoGlobalCoordinates(node_no_t nodeNoLocal)
  // global_no_t functionSpace->getNodeNoGlobalNatural(global_no_t elementNoGlobalNatural, int nodeIndex) const
  // node_no_t functionSpace->getNeighbourNodeNoLocal(node_no_t nodeNoLocal, Mesh::face_t direction)
  // node_no_t functionSpace->getNodeNo(element_no_t elementNo, int nodeIndex)
  // meshPartition->getNodeNoGlobalNatural(coordinatesGlobal)

  // run serial problem
  std::string pythonConfig = R"(

nx = 13   # number of elements in x direction
ny = 21   # number of elements in y direction

config = {
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "nElements": [nx, ny],
        "physicalExtent": [1.0, 1.0],
      }
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  typedef SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<2>,
    BasisFunction::LagrangeOfOrder<1>,
    Quadrature::Gauss<1>,
    Equation::None
  > ProblemType;
  ProblemType problem(settings);

  problem.initialize();

  LOG(INFO) << " ===================== run dof no tests ===================== ";

  typedef typename ProblemType::FunctionSpace FunctionSpace;
  const int D = FunctionSpace::dim();

  std::shared_ptr<FunctionSpace> functionSpace = problem.data().functionSpace();
  std::shared_ptr<Partition::MeshPartition<FunctionSpace>> meshPartition = functionSpace->meshPartition();

  // loop over local elements with ghosts
  for (element_no_t elementNoLocal = 0; elementNoLocal < functionSpace->nElementsLocal(); elementNoLocal++)
  {
    global_no_t elementNoGlobalNatural = meshPartition->getElementNoGlobalNatural(elementNoLocal);

    LOG(INFO) << "";
    LOG(INFO) << "elementNoLocal = " << elementNoLocal << ", elementNoGlobalNatural = " << elementNoGlobalNatural;

    // loop over local nodes, with ghosts
    for (int nodeIndex = 0; nodeIndex < functionSpace->nNodesPerElement(); nodeIndex++)
    {
      node_no_t nodeNoLocal = functionSpace->getNodeNo(elementNoLocal, nodeIndex);

      LOG(INFO) << "elementNoLocal = " << elementNoLocal << ", nodeIndex = " << nodeIndex << ", nodeNoLocal = " << nodeNoLocal;

      std::array<int,D> coordinatesGlobal = meshPartition->getNodeNoGlobalCoordinates(nodeNoLocal);

      LOG(INFO) << "coordinatesGlobal: " << coordinatesGlobal;

      global_no_t nodeNoGlobalNaturalFromCoordinates = meshPartition->getNodeNoGlobalNatural(coordinatesGlobal);
      global_no_t nodeNoGlobalNatural = functionSpace->getNodeNoGlobalNatural(elementNoGlobalNatural, nodeIndex);

      LOG(INFO) << "nodeNoGlobalNaturalFromCoordinates: " << nodeNoGlobalNaturalFromCoordinates << ", nodeNoGlobalNatural: " << nodeNoGlobalNatural;

      ASSERT_EQ(nodeNoGlobalNaturalFromCoordinates, nodeNoGlobalNatural);
    }
  }

  LOG(INFO) << " ===================== run neighbour ===================== ";
  // validate getNeighbourNodeNoLocal
  std::vector<std::pair<std::pair<Mesh::face_t,Mesh::face_t>,std::array<int,D>>> complementFaces = {
    {{Mesh::face_t::face0Plus,  Mesh::face_t::face0Minus}, {1,0}},
    {{Mesh::face_t::face0Minus, Mesh::face_t::face0Plus}, {-1,0}},
    {{Mesh::face_t::face1Plus,  Mesh::face_t::face1Minus}, {0,1}},
    {{Mesh::face_t::face1Minus, Mesh::face_t::face1Plus}, {0,-1}}
  };

  for (int i = 0; i < complementFaces.size(); i++)
  {
    Mesh::face_t face0 = complementFaces[i].first.first;
    Mesh::face_t face1 = complementFaces[i].first.second;
    std::array<int,D> expectedDifference = complementFaces[i].second;

    // loop over local nodes without ghosts
    for (node_no_t nodeNoLocal = 0; nodeNoLocal < functionSpace->nNodesLocalWithoutGhosts(); nodeNoLocal++)
    {
      node_no_t neighbourNodeNo = functionSpace->getNeighbourNodeNoLocal(nodeNoLocal, face0);
      LOG(INFO) << "nodeNoLocal = " << nodeNoLocal << ", face " << face0 << ", neighbour: " << neighbourNodeNo;

      if (neighbourNodeNo != -1)
      {
        // compare difference between global coordinates of nodes
        std::array<int,D> coordinatesGlobalNeighbour = meshPartition->getNodeNoGlobalCoordinates(neighbourNodeNo);
        std::array<int,D> coordinatesGlobalOrigin = meshPartition->getNodeNoGlobalCoordinates(nodeNoLocal);


        std::array<int,D> realDifference;
        for(int j = 0; j < D; j++)
        {
          realDifference[j] = coordinatesGlobalNeighbour[j] - coordinatesGlobalOrigin[j];
        }

        LOG(INFO) << "coordinates origin: " << coordinatesGlobalOrigin << ", neighbour: " << coordinatesGlobalNeighbour
          << ", difference: " << realDifference << ", expected: " << expectedDifference;

        for(int j = 0; j < D; j++)
        {
          EXPECT_EQ(expectedDifference[j], realDifference[j]);
        }
      }

      // if the neighbour exists and is a non-ghost dof
      if (neighbourNodeNo != -1 && neighbourNodeNo < functionSpace->nNodesLocalWithoutGhosts())
      {
        node_no_t originNode = functionSpace->getNeighbourNodeNoLocal(neighbourNodeNo, face1);

        LOG(INFO) << "  originNode: " << originNode << " face " << face1;
        ASSERT_EQ(nodeNoLocal, originNode);
      }
    }
  }
}

TEST(NumberingsTest, DofNodeNumberingDeformableQuadratic2D)
{
  // methods that are tested in this test case:
  // global_no_t meshPartition->getElementNoGlobalNatural(element_no_t elementNoLocal)
  // std::array<int,D> meshPartition_->getNodeNoGlobalCoordinates(node_no_t nodeNoLocal)
  // global_no_t functionSpace->getNodeNoGlobalNatural(global_no_t elementNoGlobalNatural, int nodeIndex) const
  // node_no_t functionSpace->getNeighbourNodeNoLocal(node_no_t nodeNoLocal, Mesh::face_t direction)
  // node_no_t functionSpace->getNodeNo(element_no_t elementNo, int nodeIndex)
  // meshPartition->getNodeNoGlobalNatural(coordinatesGlobal)

  // run serial problem
  std::string pythonConfig = R"(

nx = 13   # number of elements in x direction
ny = 21   # number of elements in y direction

config = {
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "nElements": [nx, ny],
        "physicalExtent": [1.0, 1.0],
      }
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  typedef SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<2>,
    BasisFunction::LagrangeOfOrder<2>,
    Quadrature::Gauss<1>,
    Equation::None
  > ProblemType;
  ProblemType problem(settings);

  problem.initialize();

  LOG(INFO) << " ===================== run dof no tests ===================== ";

  typedef typename ProblemType::FunctionSpace FunctionSpace;
  const int D = FunctionSpace::dim();

  std::shared_ptr<FunctionSpace> functionSpace = problem.data().functionSpace();
  std::shared_ptr<Partition::MeshPartition<FunctionSpace>> meshPartition = functionSpace->meshPartition();

  // loop over local elements with ghosts
  for (element_no_t elementNoLocal = 0; elementNoLocal < functionSpace->nElementsLocal(); elementNoLocal++)
  {
    global_no_t elementNoGlobalNatural = meshPartition->getElementNoGlobalNatural(elementNoLocal);

    LOG(INFO) << "";
    LOG(INFO) << "elementNoLocal = " << elementNoLocal << ", elementNoGlobalNatural = " << elementNoGlobalNatural;

    // loop over local nodes, with ghosts
    for (int nodeIndex = 0; nodeIndex < functionSpace->nNodesPerElement(); nodeIndex++)
    {
      node_no_t nodeNoLocal = functionSpace->getNodeNo(elementNoLocal, nodeIndex);

      LOG(INFO) << "elementNoLocal = " << elementNoLocal << ", nodeIndex = " << nodeIndex << ", nodeNoLocal = " << nodeNoLocal;

      std::array<int,D> coordinatesGlobal = meshPartition->getNodeNoGlobalCoordinates(nodeNoLocal);

      LOG(INFO) << "coordinatesGlobal: " << coordinatesGlobal;

      global_no_t nodeNoGlobalNaturalFromCoordinates = meshPartition->getNodeNoGlobalNatural(coordinatesGlobal);
      global_no_t nodeNoGlobalNatural = functionSpace->getNodeNoGlobalNatural(elementNoGlobalNatural, nodeIndex);

      LOG(INFO) << "nodeNoGlobalNaturalFromCoordinates: " << nodeNoGlobalNaturalFromCoordinates << ", nodeNoGlobalNatural: " << nodeNoGlobalNatural;

      ASSERT_EQ(nodeNoGlobalNaturalFromCoordinates, nodeNoGlobalNatural);
    }
  }

  LOG(INFO) << " ===================== run neighbour ===================== ";
  // validate getNeighbourNodeNoLocal
  std::vector<std::pair<std::pair<Mesh::face_t,Mesh::face_t>,std::array<int,D>>> complementFaces = {
    {{Mesh::face_t::face0Plus,  Mesh::face_t::face0Minus}, {1,0}},
    {{Mesh::face_t::face0Minus, Mesh::face_t::face0Plus}, {-1,0}},
    {{Mesh::face_t::face1Plus,  Mesh::face_t::face1Minus}, {0,1}},
    {{Mesh::face_t::face1Minus, Mesh::face_t::face1Plus}, {0,-1}}
  };

  for (int i = 0; i < complementFaces.size(); i++)
  {
    Mesh::face_t face0 = complementFaces[i].first.first;
    Mesh::face_t face1 = complementFaces[i].first.second;
    std::array<int,D> expectedDifference = complementFaces[i].second;

    // loop over local nodes without ghosts
    for (node_no_t nodeNoLocal = 0; nodeNoLocal < functionSpace->nNodesLocalWithoutGhosts(); nodeNoLocal++)
    {
      node_no_t neighbourNodeNo = functionSpace->getNeighbourNodeNoLocal(nodeNoLocal, face0);
      LOG(INFO) << "nodeNoLocal = " << nodeNoLocal << ", face " << face0 << ", neighbour: " << neighbourNodeNo;

      if (neighbourNodeNo != -1)
      {
        // compare difference between global coordinates of nodes
        std::array<int,D> coordinatesGlobalNeighbour = meshPartition->getNodeNoGlobalCoordinates(neighbourNodeNo);
        std::array<int,D> coordinatesGlobalOrigin = meshPartition->getNodeNoGlobalCoordinates(nodeNoLocal);


        std::array<int,D> realDifference;
        for(int j = 0; j < D; j++)
        {
          realDifference[j] = coordinatesGlobalNeighbour[j] - coordinatesGlobalOrigin[j];
        }

        LOG(INFO) << "coordinates origin: " << coordinatesGlobalOrigin << ", neighbour: " << coordinatesGlobalNeighbour
          << ", difference: " << realDifference << ", expected: " << expectedDifference;

        for(int j = 0; j < D; j++)
        {
          EXPECT_EQ(expectedDifference[j], realDifference[j]);
        }
      }

      // if the neighbour exists and is a non-ghost dof
      if (neighbourNodeNo != -1 && neighbourNodeNo < functionSpace->nNodesLocalWithoutGhosts())
      {
        node_no_t originNode = functionSpace->getNeighbourNodeNoLocal(neighbourNodeNo, face1);

        LOG(INFO) << "  originNode: " << originNode << " face " << face1;
        ASSERT_EQ(nodeNoLocal, originNode);
      }
    }
  }
}

TEST(NumberingsTest, DofNodeNumberingDeformableHermite2D)
{
  // methods that are tested in this test case:
  // global_no_t meshPartition->getElementNoGlobalNatural(element_no_t elementNoLocal)
  // std::array<int,D> meshPartition_->getNodeNoGlobalCoordinates(node_no_t nodeNoLocal)
  // global_no_t functionSpace->getNodeNoGlobalNatural(global_no_t elementNoGlobalNatural, int nodeIndex) const
  // node_no_t functionSpace->getNeighbourNodeNoLocal(node_no_t nodeNoLocal, Mesh::face_t direction)
  // node_no_t functionSpace->getNodeNo(element_no_t elementNo, int nodeIndex)
  // meshPartition->getNodeNoGlobalNatural(coordinatesGlobal)

  // run serial problem
  std::string pythonConfig = R"(

nx = 13   # number of elements in x direction
ny = 21   # number of elements in y direction

config = {
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "nElements": [nx, ny],
        "physicalExtent": [1.0, 1.0],
      }
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  typedef SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<2>,
    BasisFunction::Hermite,
    Quadrature::Gauss<1>,
    Equation::None
  > ProblemType;
  ProblemType problem(settings);

  problem.initialize();

  LOG(INFO) << " ===================== run dof no tests ===================== ";

  typedef typename ProblemType::FunctionSpace FunctionSpace;
  const int D = FunctionSpace::dim();

  std::shared_ptr<FunctionSpace> functionSpace = problem.data().functionSpace();
  std::shared_ptr<Partition::MeshPartition<FunctionSpace>> meshPartition = functionSpace->meshPartition();

  // loop over local elements with ghosts
  for (element_no_t elementNoLocal = 0; elementNoLocal < functionSpace->nElementsLocal(); elementNoLocal++)
  {
    global_no_t elementNoGlobalNatural = meshPartition->getElementNoGlobalNatural(elementNoLocal);

    LOG(INFO) << "";
    LOG(INFO) << "elementNoLocal = " << elementNoLocal << ", elementNoGlobalNatural = " << elementNoGlobalNatural;

    // loop over local nodes, with ghosts
    for (int nodeIndex = 0; nodeIndex < functionSpace->nNodesPerElement(); nodeIndex++)
    {
      node_no_t nodeNoLocal = functionSpace->getNodeNo(elementNoLocal, nodeIndex);

      LOG(INFO) << "elementNoLocal = " << elementNoLocal << ", nodeIndex = " << nodeIndex << ", nodeNoLocal = " << nodeNoLocal;

      std::array<int,D> coordinatesGlobal = meshPartition->getNodeNoGlobalCoordinates(nodeNoLocal);

      LOG(INFO) << "coordinatesGlobal: " << coordinatesGlobal;

      global_no_t nodeNoGlobalNaturalFromCoordinates = meshPartition->getNodeNoGlobalNatural(coordinatesGlobal);
      global_no_t nodeNoGlobalNatural = functionSpace->getNodeNoGlobalNatural(elementNoGlobalNatural, nodeIndex);

      LOG(INFO) << "nodeNoGlobalNaturalFromCoordinates: " << nodeNoGlobalNaturalFromCoordinates << ", nodeNoGlobalNatural: " << nodeNoGlobalNatural;

      ASSERT_EQ(nodeNoGlobalNaturalFromCoordinates, nodeNoGlobalNatural);
    }
  }

  LOG(INFO) << " ===================== run neighbour ===================== ";
  // validate getNeighbourNodeNoLocal
  std::vector<std::pair<Mesh::face_t,Mesh::face_t>> complementFaces = {
    {Mesh::face_t::face0Plus,  Mesh::face_t::face0Minus},
    {Mesh::face_t::face0Minus, Mesh::face_t::face0Plus},
    {Mesh::face_t::face1Plus,  Mesh::face_t::face1Minus},
    {Mesh::face_t::face1Minus, Mesh::face_t::face1Plus},
  };

  for (int i = 0; i < complementFaces.size(); i++)
  {
    Mesh::face_t face0 = complementFaces[i].first;
    Mesh::face_t face1 = complementFaces[i].second;

    // loop over local nodes without ghosts
    for (node_no_t nodeNoLocal = 0; nodeNoLocal < functionSpace->nNodesLocalWithoutGhosts(); nodeNoLocal++)
    {
      node_no_t neighbourNodeNo = functionSpace->getNeighbourNodeNoLocal(nodeNoLocal, face0);
      LOG(INFO) << "nodeNoLocal = " << nodeNoLocal << ", face " << face0 << ", neighbour: " << neighbourNodeNo;

      // if the neighbour exists and is a non-ghost dof
      if (neighbourNodeNo != -1 && neighbourNodeNo < functionSpace->nNodesLocalWithoutGhosts())
      {
        node_no_t originNode = functionSpace->getNeighbourNodeNoLocal(neighbourNodeNo, face1);

        LOG(INFO) << "  originNode: " << originNode << " face " << face1;
        ASSERT_EQ(nodeNoLocal, originNode);
      }
    }
  }
}


// ---------- 3D ---------------
TEST(NumberingsTest, DofNodeNumberingDeformableLinear3D)
{
  // methods that are tested in this test case:
  // global_no_t meshPartition->getElementNoGlobalNatural(element_no_t elementNoLocal)
  // std::array<int,D> meshPartition_->getNodeNoGlobalCoordinates(node_no_t nodeNoLocal)
  // global_no_t functionSpace->getNodeNoGlobalNatural(global_no_t elementNoGlobalNatural, int nodeIndex) const
  // node_no_t functionSpace->getNeighbourNodeNoLocal(node_no_t nodeNoLocal, Mesh::face_t direction)
  // node_no_t functionSpace->getNodeNo(element_no_t elementNo, int nodeIndex)
  // meshPartition->getNodeNoGlobalNatural(coordinatesGlobal)

  // run serial problem
  std::string pythonConfig = R"(

nx = 11   # number of elements in x direction
ny = 17   # number of elements in y direction
nz = 6   # number of elements in z direction

config = {
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "nElements": [nx, ny, nz],
        "physicalExtent": [1.0, 1.0, 1.0],
      }
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  typedef SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<1>,
    Quadrature::Gauss<1>,
    Equation::None
  > ProblemType;
  ProblemType problem(settings);

  problem.initialize();

  LOG(INFO) << " ===================== run dof no tests ===================== ";

  typedef typename ProblemType::FunctionSpace FunctionSpace;
  const int D = FunctionSpace::dim();

  std::shared_ptr<FunctionSpace> functionSpace = problem.data().functionSpace();
  std::shared_ptr<Partition::MeshPartition<FunctionSpace>> meshPartition = functionSpace->meshPartition();

  // loop over local elements with ghosts
  for (element_no_t elementNoLocal = 0; elementNoLocal < functionSpace->nElementsLocal(); elementNoLocal++)
  {
    global_no_t elementNoGlobalNatural = meshPartition->getElementNoGlobalNatural(elementNoLocal);

    LOG(INFO) << "";
    LOG(INFO) << "elementNoLocal = " << elementNoLocal << ", elementNoGlobalNatural = " << elementNoGlobalNatural;

    // loop over local nodes, with ghosts
    for (int nodeIndex = 0; nodeIndex < functionSpace->nNodesPerElement(); nodeIndex++)
    {
      node_no_t nodeNoLocal = functionSpace->getNodeNo(elementNoLocal, nodeIndex);

      LOG(INFO) << "elementNoLocal = " << elementNoLocal << ", nodeIndex = " << nodeIndex << ", nodeNoLocal = " << nodeNoLocal;

      std::array<int,D> coordinatesGlobal = meshPartition->getNodeNoGlobalCoordinates(nodeNoLocal);

      LOG(INFO) << "coordinatesGlobal: " << coordinatesGlobal;

      global_no_t nodeNoGlobalNaturalFromCoordinates = meshPartition->getNodeNoGlobalNatural(coordinatesGlobal);
      global_no_t nodeNoGlobalNatural = functionSpace->getNodeNoGlobalNatural(elementNoGlobalNatural, nodeIndex);

      LOG(INFO) << "nodeNoGlobalNaturalFromCoordinates: " << nodeNoGlobalNaturalFromCoordinates << ", nodeNoGlobalNatural: " << nodeNoGlobalNatural;

      ASSERT_EQ(nodeNoGlobalNaturalFromCoordinates, nodeNoGlobalNatural);
    }
  }

  LOG(INFO) << " ===================== run neighbour ===================== ";
  // validate getNeighbourNodeNoLocal
  std::vector<std::pair<std::pair<Mesh::face_t,Mesh::face_t>,std::array<int,D>>> complementFaces = {
    {{Mesh::face_t::face0Plus,  Mesh::face_t::face0Minus}, {1,0,0}},
    {{Mesh::face_t::face0Minus, Mesh::face_t::face0Plus}, {-1,0,0}},
    {{Mesh::face_t::face1Plus,  Mesh::face_t::face1Minus}, {0,1,0}},
    {{Mesh::face_t::face1Minus, Mesh::face_t::face1Plus}, {0,-1,0}},
    {{Mesh::face_t::face2Plus,  Mesh::face_t::face2Minus}, {0,0,1}},
    {{Mesh::face_t::face2Minus, Mesh::face_t::face2Plus}, {0,0,-1}}
  };

  for (int i = 0; i < complementFaces.size(); i++)
  {
    Mesh::face_t face0 = complementFaces[i].first.first;
    Mesh::face_t face1 = complementFaces[i].first.second;
    std::array<int,D> expectedDifference = complementFaces[i].second;

    // loop over local nodes without ghosts
    for (node_no_t nodeNoLocal = 0; nodeNoLocal < functionSpace->nNodesLocalWithoutGhosts(); nodeNoLocal++)
    {
      node_no_t neighbourNodeNo = functionSpace->getNeighbourNodeNoLocal(nodeNoLocal, face0);
      LOG(INFO) << "nodeNoLocal = " << nodeNoLocal << ", face " << face0 << ", neighbour: " << neighbourNodeNo;

      if (neighbourNodeNo != -1)
      {
        // compare difference between global coordinates of nodes
        std::array<int,D> coordinatesGlobalNeighbour = meshPartition->getNodeNoGlobalCoordinates(neighbourNodeNo);
        std::array<int,D> coordinatesGlobalOrigin = meshPartition->getNodeNoGlobalCoordinates(nodeNoLocal);


        std::array<int,D> realDifference;
        for(int j = 0; j < D; j++)
        {
          realDifference[j] = coordinatesGlobalNeighbour[j] - coordinatesGlobalOrigin[j];
        }

        LOG(INFO) << "coordinates origin: " << coordinatesGlobalOrigin << ", neighbour: " << coordinatesGlobalNeighbour
          << ", difference: " << realDifference << ", expected: " << expectedDifference;

        for(int j = 0; j < D; j++)
        {
          EXPECT_EQ(expectedDifference[j], realDifference[j]);
        }
      }

      // if the neighbour exists and is a non-ghost dof
      if (neighbourNodeNo != -1 && neighbourNodeNo < functionSpace->nNodesLocalWithoutGhosts())
      {
        node_no_t originNode = functionSpace->getNeighbourNodeNoLocal(neighbourNodeNo, face1);

        LOG(INFO) << "  originNode: " << originNode << " face " << face1;
        ASSERT_EQ(nodeNoLocal, originNode);
      }
    }
  }
}

TEST(NumberingsTest, DofNodeNumberingDeformableQuadratic3D)
{
  // methods that are tested in this test case:
  // global_no_t meshPartition->getElementNoGlobalNatural(element_no_t elementNoLocal)
  // std::array<int,D> meshPartition_->getNodeNoGlobalCoordinates(node_no_t nodeNoLocal)
  // global_no_t functionSpace->getNodeNoGlobalNatural(global_no_t elementNoGlobalNatural, int nodeIndex) const
  // node_no_t functionSpace->getNeighbourNodeNoLocal(node_no_t nodeNoLocal, Mesh::face_t direction)
  // node_no_t functionSpace->getNodeNo(element_no_t elementNo, int nodeIndex)
  // meshPartition->getNodeNoGlobalNatural(coordinatesGlobal)

  // run serial problem
  std::string pythonConfig = R"(

nx = 11   # number of elements in x direction
ny = 17   # number of elements in y direction
nz = 6   # number of elements in z direction

config = {
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "nElements": [nx, ny, nz],
        "physicalExtent": [1.0, 1.0, 1.0],
      }
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  typedef SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::LagrangeOfOrder<2>,
    Quadrature::Gauss<1>,
    Equation::None
  > ProblemType;
  ProblemType problem(settings);

  problem.initialize();

  LOG(INFO) << " ===================== run dof no tests ===================== ";

  typedef typename ProblemType::FunctionSpace FunctionSpace;
  const int D = FunctionSpace::dim();

  std::shared_ptr<FunctionSpace> functionSpace = problem.data().functionSpace();
  std::shared_ptr<Partition::MeshPartition<FunctionSpace>> meshPartition = functionSpace->meshPartition();

  // loop over local elements with ghosts
  for (element_no_t elementNoLocal = 0; elementNoLocal < functionSpace->nElementsLocal(); elementNoLocal++)
  {
    global_no_t elementNoGlobalNatural = meshPartition->getElementNoGlobalNatural(elementNoLocal);

    LOG(INFO) << "";
    LOG(INFO) << "elementNoLocal = " << elementNoLocal << ", elementNoGlobalNatural = " << elementNoGlobalNatural;

    // loop over local nodes, with ghosts
    for (int nodeIndex = 0; nodeIndex < functionSpace->nNodesPerElement(); nodeIndex++)
    {
      node_no_t nodeNoLocal = functionSpace->getNodeNo(elementNoLocal, nodeIndex);

      LOG(INFO) << "elementNoLocal = " << elementNoLocal << ", nodeIndex = " << nodeIndex << ", nodeNoLocal = " << nodeNoLocal;

      std::array<int,D> coordinatesGlobal = meshPartition->getNodeNoGlobalCoordinates(nodeNoLocal);

      LOG(INFO) << "coordinatesGlobal: " << coordinatesGlobal;

      global_no_t nodeNoGlobalNaturalFromCoordinates = meshPartition->getNodeNoGlobalNatural(coordinatesGlobal);
      global_no_t nodeNoGlobalNatural = functionSpace->getNodeNoGlobalNatural(elementNoGlobalNatural, nodeIndex);

      LOG(INFO) << "nodeNoGlobalNaturalFromCoordinates: " << nodeNoGlobalNaturalFromCoordinates << ", nodeNoGlobalNatural: " << nodeNoGlobalNatural;

      ASSERT_EQ(nodeNoGlobalNaturalFromCoordinates, nodeNoGlobalNatural);
    }
  }

  LOG(INFO) << " ===================== run neighbour ===================== ";
  // validate getNeighbourNodeNoLocal
  std::vector<std::pair<std::pair<Mesh::face_t,Mesh::face_t>,std::array<int,D>>> complementFaces = {
    {{Mesh::face_t::face0Plus,  Mesh::face_t::face0Minus}, {1,0,0}},
    {{Mesh::face_t::face0Minus, Mesh::face_t::face0Plus}, {-1,0,0}},
    {{Mesh::face_t::face1Plus,  Mesh::face_t::face1Minus}, {0,1,0}},
    {{Mesh::face_t::face1Minus, Mesh::face_t::face1Plus}, {0,-1,0}},
    {{Mesh::face_t::face2Plus,  Mesh::face_t::face2Minus}, {0,0,1}},
    {{Mesh::face_t::face2Minus, Mesh::face_t::face2Plus}, {0,0,-1}}
  };

  for (int i = 0; i < complementFaces.size(); i++)
  {
    Mesh::face_t face0 = complementFaces[i].first.first;
    Mesh::face_t face1 = complementFaces[i].first.second;
    std::array<int,D> expectedDifference = complementFaces[i].second;

    // loop over local nodes without ghosts
    for (node_no_t nodeNoLocal = 0; nodeNoLocal < functionSpace->nNodesLocalWithoutGhosts(); nodeNoLocal++)
    {
      node_no_t neighbourNodeNo = functionSpace->getNeighbourNodeNoLocal(nodeNoLocal, face0);
      LOG(INFO) << "nodeNoLocal = " << nodeNoLocal << ", face " << face0 << ", neighbour: " << neighbourNodeNo;

      if (neighbourNodeNo != -1)
      {
        // compare difference between global coordinates of nodes
        std::array<int,D> coordinatesGlobalNeighbour = meshPartition->getNodeNoGlobalCoordinates(neighbourNodeNo);
        std::array<int,D> coordinatesGlobalOrigin = meshPartition->getNodeNoGlobalCoordinates(nodeNoLocal);


        std::array<int,D> realDifference;
        for(int j = 0; j < D; j++)
        {
          realDifference[j] = coordinatesGlobalNeighbour[j] - coordinatesGlobalOrigin[j];
        }

        LOG(INFO) << "coordinates origin: " << coordinatesGlobalOrigin << ", neighbour: " << coordinatesGlobalNeighbour
          << ", difference: " << realDifference << ", expected: " << expectedDifference;

        for(int j = 0; j < D; j++)
        {
          EXPECT_EQ(expectedDifference[j], realDifference[j]);
        }
      }

      // if the neighbour exists and is a non-ghost dof
      if (neighbourNodeNo != -1 && neighbourNodeNo < functionSpace->nNodesLocalWithoutGhosts())
      {
        node_no_t originNode = functionSpace->getNeighbourNodeNoLocal(neighbourNodeNo, face1);

        LOG(INFO) << "  originNode: " << originNode << " face " << face1;
        ASSERT_EQ(nodeNoLocal, originNode);
      }
    }
  }
}

TEST(NumberingsTest, DofNodeNumberingDeformableHermite3D)
{
  // methods that are tested in this test case:
  // global_no_t meshPartition->getElementNoGlobalNatural(element_no_t elementNoLocal)
  // std::array<int,D> meshPartition_->getNodeNoGlobalCoordinates(node_no_t nodeNoLocal)
  // global_no_t functionSpace->getNodeNoGlobalNatural(global_no_t elementNoGlobalNatural, int nodeIndex) const
  // node_no_t functionSpace->getNeighbourNodeNoLocal(node_no_t nodeNoLocal, Mesh::face_t direction)
  // node_no_t functionSpace->getNodeNo(element_no_t elementNo, int nodeIndex)
  // meshPartition->getNodeNoGlobalNatural(coordinatesGlobal)

  // run serial problem
  std::string pythonConfig = R"(

nx = 11   # number of elements in x direction
ny = 17   # number of elements in y direction
nz = 6   # number of elements in z direction

config = {
      "FiniteElementMethod": {
        "inputMeshIsGlobal": True,
        "nElements": [nx, ny, nz],
        "physicalExtent": [1.0, 1.0, 1.0],
      }
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  typedef SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunction::Hermite,
    Quadrature::Gauss<1>,
    Equation::None
  > ProblemType;
  ProblemType problem(settings);

  problem.initialize();

  LOG(INFO) << " ===================== run dof no tests ===================== ";

  typedef typename ProblemType::FunctionSpace FunctionSpace;
  const int D = FunctionSpace::dim();

  std::shared_ptr<FunctionSpace> functionSpace = problem.data().functionSpace();
  std::shared_ptr<Partition::MeshPartition<FunctionSpace>> meshPartition = functionSpace->meshPartition();

  // loop over local elements with ghosts
  for (element_no_t elementNoLocal = 0; elementNoLocal < functionSpace->nElementsLocal(); elementNoLocal++)
  {
    global_no_t elementNoGlobalNatural = meshPartition->getElementNoGlobalNatural(elementNoLocal);

    LOG(INFO) << "";
    LOG(INFO) << "elementNoLocal = " << elementNoLocal << ", elementNoGlobalNatural = " << elementNoGlobalNatural;

    // loop over local nodes, with ghosts
    for (int nodeIndex = 0; nodeIndex < functionSpace->nNodesPerElement(); nodeIndex++)
    {
      node_no_t nodeNoLocal = functionSpace->getNodeNo(elementNoLocal, nodeIndex);

      LOG(INFO) << "elementNoLocal = " << elementNoLocal << ", nodeIndex = " << nodeIndex << ", nodeNoLocal = " << nodeNoLocal;

      std::array<int,D> coordinatesGlobal = meshPartition->getNodeNoGlobalCoordinates(nodeNoLocal);

      LOG(INFO) << "coordinatesGlobal: " << coordinatesGlobal;

      global_no_t nodeNoGlobalNaturalFromCoordinates = meshPartition->getNodeNoGlobalNatural(coordinatesGlobal);
      global_no_t nodeNoGlobalNatural = functionSpace->getNodeNoGlobalNatural(elementNoGlobalNatural, nodeIndex);

      LOG(INFO) << "nodeNoGlobalNaturalFromCoordinates: " << nodeNoGlobalNaturalFromCoordinates << ", nodeNoGlobalNatural: " << nodeNoGlobalNatural;

      ASSERT_EQ(nodeNoGlobalNaturalFromCoordinates, nodeNoGlobalNatural);
    }
  }

  LOG(INFO) << " ===================== run neighbour ===================== ";
  // validate getNeighbourNodeNoLocal
  std::vector<std::pair<std::pair<Mesh::face_t,Mesh::face_t>,std::array<int,D>>> complementFaces = {
    {{Mesh::face_t::face0Plus,  Mesh::face_t::face0Minus}, {1,0,0}},
    {{Mesh::face_t::face0Minus, Mesh::face_t::face0Plus}, {-1,0,0}},
    {{Mesh::face_t::face1Plus,  Mesh::face_t::face1Minus}, {0,1,0}},
    {{Mesh::face_t::face1Minus, Mesh::face_t::face1Plus}, {0,-1,0}},
    {{Mesh::face_t::face2Plus,  Mesh::face_t::face2Minus}, {0,0,1}},
    {{Mesh::face_t::face2Minus, Mesh::face_t::face2Plus}, {0,0,-1}}
  };

  for (int i = 0; i < complementFaces.size(); i++)
  {
    Mesh::face_t face0 = complementFaces[i].first.first;
    Mesh::face_t face1 = complementFaces[i].first.second;
    std::array<int,D> expectedDifference = complementFaces[i].second;

    // loop over local nodes without ghosts
    for (node_no_t nodeNoLocal = 0; nodeNoLocal < functionSpace->nNodesLocalWithoutGhosts(); nodeNoLocal++)
    {
      node_no_t neighbourNodeNo = functionSpace->getNeighbourNodeNoLocal(nodeNoLocal, face0);
      LOG(INFO) << "nodeNoLocal = " << nodeNoLocal << ", face " << face0 << ", neighbour: " << neighbourNodeNo;

      if (neighbourNodeNo != -1)
      {
        // compare difference between global coordinates of nodes
        std::array<int,D> coordinatesGlobalNeighbour = meshPartition->getNodeNoGlobalCoordinates(neighbourNodeNo);
        std::array<int,D> coordinatesGlobalOrigin = meshPartition->getNodeNoGlobalCoordinates(nodeNoLocal);


        std::array<int,D> realDifference;
        for(int j = 0; j < D; j++)
        {
          realDifference[j] = coordinatesGlobalNeighbour[j] - coordinatesGlobalOrigin[j];
        }

        LOG(INFO) << "coordinates origin: " << coordinatesGlobalOrigin << ", neighbour: " << coordinatesGlobalNeighbour
          << ", difference: " << realDifference << ", expected: " << expectedDifference;

        for(int j = 0; j < D; j++)
        {
          EXPECT_EQ(expectedDifference[j], realDifference[j]);
        }
      }

      // if the neighbour exists and is a non-ghost dof
      if (neighbourNodeNo != -1 && neighbourNodeNo < functionSpace->nNodesLocalWithoutGhosts())
      {
        node_no_t originNode = functionSpace->getNeighbourNodeNoLocal(neighbourNodeNo, face1);

        LOG(INFO) << "  originNode: " << originNode << " face " << face1;
        ASSERT_EQ(nodeNoLocal, originNode);
      }
    }
  }
}
