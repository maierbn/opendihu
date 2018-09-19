#include <Python.h>  // this has to be the first included header

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cassert>
#include <cmath>

#include "gtest/gtest.h"
#include "opendihu.h"
#include "arg.h"
#include "stiffness_matrix_tester.h"
#include "node_positions_tester.h"
#include "../utility.h"
#include "mesh/face_t.h"
#include "function_space/function_space.h"

namespace SpatialDiscretization
{
  
  
TEST(FaceTest, parseFaceString)
{ 
  ASSERT_EQ(Mesh::parseFace("0-"), Mesh::face_t::face0Minus);
  ASSERT_EQ(Mesh::parseFace("0+"), Mesh::face_t::face0Plus);
  ASSERT_EQ(Mesh::parseFace("1-"), Mesh::face_t::face1Minus);
  ASSERT_EQ(Mesh::parseFace("1+"), Mesh::face_t::face1Plus);
  ASSERT_EQ(Mesh::parseFace("2-"), Mesh::face_t::face2Minus);
  ASSERT_EQ(Mesh::parseFace("2+"), Mesh::face_t::face2Plus);
}
  
TEST(FaceTest, faceDofsLinearLagrange1D)
{ 
  typedef Mesh::UnstructuredDeformableOfDimension<1> MeshType;
  typedef BasisFunction::LagrangeOfOrder<1> BasisFunctionType;
  typedef FunctionSpace::FunctionSpace<MeshType, BasisFunctionType> FunctionSpaceType;
  
  std::array<dof_no_t, 1> dofIndices1;
  std::array<dof_no_t, 1> dofIndices1Reference = {0};
  FunctionSpaceType::getFaceDofs(Mesh::face_t::face0Minus, dofIndices1);  // 0-
  ASSERT_EQ(dofIndices1, dofIndices1Reference);
  
  std::array<dof_no_t, 1> dofIndices2;
  std::array<dof_no_t, 1> dofIndices2Reference = {1};
  FunctionSpaceType::getFaceDofs(Mesh::face_t::face0Plus, dofIndices2);  // 0+
  ASSERT_EQ(dofIndices2, dofIndices2Reference);
}
  
TEST(FaceTest, faceDofsQuadraticLagrange1D)
{ 
  typedef Mesh::UnstructuredDeformableOfDimension<1> MeshType;
  typedef BasisFunction::LagrangeOfOrder<2> BasisFunctionType;
  typedef FunctionSpace::FunctionSpace<MeshType, BasisFunctionType> FunctionSpaceType;
  
  std::array<dof_no_t, 1> dofIndices1;
  std::array<dof_no_t, 1> dofIndices1Reference = {0};
  FunctionSpaceType::getFaceDofs(Mesh::face_t::face0Minus, dofIndices1);  // 0-
  ASSERT_EQ(dofIndices1, dofIndices1Reference);
  
  std::array<dof_no_t, 1> dofIndices2;
  std::array<dof_no_t, 1> dofIndices2Reference = {2};
  FunctionSpaceType::getFaceDofs(Mesh::face_t::face0Plus, dofIndices2);  // 0+
  ASSERT_EQ(dofIndices2, dofIndices2Reference);
}

TEST(FaceTest, faceDofsLinearLagrange2D)
{ 
  typedef Mesh::UnstructuredDeformableOfDimension<2> MeshType;
  typedef BasisFunction::LagrangeOfOrder<1> BasisFunctionType;
  typedef FunctionSpace::FunctionSpace<MeshType, BasisFunctionType> FunctionSpaceType;
  
  std::array<dof_no_t, 2> dofIndices1;
  std::array<dof_no_t, 2> dofIndices1Reference = {0, 2};
  FunctionSpaceType::getFaceDofs(Mesh::face_t::face0Minus, dofIndices1);  // 0-
  ASSERT_EQ(dofIndices1, dofIndices1Reference);
  
  std::array<dof_no_t, 2> dofIndices2;
  std::array<dof_no_t, 2> dofIndices2Reference = {1, 3};
  FunctionSpaceType::getFaceDofs(Mesh::face_t::face0Plus, dofIndices2);  // 0+
  ASSERT_EQ(dofIndices2, dofIndices2Reference);
  
  std::array<dof_no_t, 2> dofIndices3;
  std::array<dof_no_t, 2> dofIndices3Reference = {0, 1};
  FunctionSpaceType::getFaceDofs(Mesh::face_t::face1Minus, dofIndices3);  // 1-
  ASSERT_EQ(dofIndices3, dofIndices3Reference);
  
  std::array<dof_no_t, 2> dofIndices4;
  std::array<dof_no_t, 2> dofIndices4Reference = {2, 3};
  FunctionSpaceType::getFaceDofs(Mesh::face_t::face1Plus, dofIndices4);  // 1+
  ASSERT_EQ(dofIndices4, dofIndices4Reference);
}

TEST(FaceTest, faceDofsQuadraticLagrange2D)
{ 
  typedef Mesh::UnstructuredDeformableOfDimension<2> MeshType;
  typedef BasisFunction::LagrangeOfOrder<2> BasisFunctionType;
  typedef FunctionSpace::FunctionSpace<MeshType, BasisFunctionType> FunctionSpaceType;
  
  std::array<dof_no_t, 3> dofIndices1;
  std::array<dof_no_t, 3> dofIndices1Reference = {0, 3, 6};
  FunctionSpaceType::getFaceDofs(Mesh::face_t::face0Minus, dofIndices1);  // 0-
  ASSERT_EQ(dofIndices1, dofIndices1Reference);
  
  std::array<dof_no_t, 3> dofIndices2;
  std::array<dof_no_t, 3> dofIndices2Reference = {2, 5, 8};
  FunctionSpaceType::getFaceDofs(Mesh::face_t::face0Plus, dofIndices2);  // 0+
  ASSERT_EQ(dofIndices2, dofIndices2Reference);
  
  std::array<dof_no_t, 3> dofIndices3;
  std::array<dof_no_t, 3> dofIndices3Reference = {0, 1, 2};
  FunctionSpaceType::getFaceDofs(Mesh::face_t::face1Minus, dofIndices3);  // 1-
  ASSERT_EQ(dofIndices3, dofIndices3Reference);
  
  std::array<dof_no_t, 3> dofIndices4;
  std::array<dof_no_t, 3> dofIndices4Reference = {6, 7, 8};
  FunctionSpaceType::getFaceDofs(Mesh::face_t::face1Plus, dofIndices4);  // 1+
  ASSERT_EQ(dofIndices4, dofIndices4Reference);
}

TEST(FaceTest, faceDofsLinearLagrange3D)
{ 
  typedef Mesh::UnstructuredDeformableOfDimension<3> MeshType;
  typedef BasisFunction::LagrangeOfOrder<1> BasisFunctionType;
  typedef FunctionSpace::FunctionSpace<MeshType, BasisFunctionType> FunctionSpaceType;
  
  /*
  std::string pythonConfig = R"(
# 3D
  
config = {
  "FiniteElementMethod" : {
    "dirichletBoundaryConditions": bc,
  }
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  FiniteElementMethod<
    MeshType,
    BasisFunction::LagrangeOfOrder<1>,
    Quadrature::Gauss<2>,
    Equation::None
  > problem(settings);
  
  
  std::shared_ptr<MeshType> functionSpace = std::static_ptr_cast<MeshType>(problem->functionSpace());*/
  
  std::array<dof_no_t, 4> dofIndices1;
  std::array<dof_no_t, 4> dofIndices1Reference = {0, 2, 4, 6};
  FunctionSpaceType::getFaceDofs(Mesh::face_t::face0Minus, dofIndices1);  // 0-
  ASSERT_EQ(dofIndices1, dofIndices1Reference);
  
  std::array<dof_no_t, 4> dofIndices2;
  std::array<dof_no_t, 4> dofIndices2Reference = {1, 3, 5, 7};
  FunctionSpaceType::getFaceDofs(Mesh::face_t::face0Plus, dofIndices2);  // 0+
  ASSERT_EQ(dofIndices2, dofIndices2Reference);
  
  std::array<dof_no_t, 4> dofIndices3;
  std::array<dof_no_t, 4> dofIndices3Reference = {0, 1, 4, 5};
  FunctionSpaceType::getFaceDofs(Mesh::face_t::face1Minus, dofIndices3);  // 1-
  ASSERT_EQ(dofIndices3, dofIndices3Reference);
  
  std::array<dof_no_t, 4> dofIndices4;
  std::array<dof_no_t, 4> dofIndices4Reference = {2, 3, 6, 7};
  FunctionSpaceType::getFaceDofs(Mesh::face_t::face1Plus, dofIndices4);  // 1+
  ASSERT_EQ(dofIndices4, dofIndices4Reference);
  
  std::array<dof_no_t, 4> dofIndices5;
  std::array<dof_no_t, 4> dofIndices5Reference = {0, 1, 2, 3};
  FunctionSpaceType::getFaceDofs(Mesh::face_t::face2Minus, dofIndices5);  // 2-
  ASSERT_EQ(dofIndices5, dofIndices5Reference);
  
  std::array<dof_no_t, 4> dofIndices6;
  std::array<dof_no_t, 4> dofIndices6Reference = {4, 5, 6, 7};
  FunctionSpaceType::getFaceDofs(Mesh::face_t::face2Plus, dofIndices6);  // 2+
  ASSERT_EQ(dofIndices6, dofIndices6Reference);
  
  
}

TEST(FaceTest, faceDofsQuadraticLagrange3D)
{ 
  typedef Mesh::UnstructuredDeformableOfDimension<3> MeshType;
  typedef BasisFunction::LagrangeOfOrder<2> BasisFunctionType;
  typedef FunctionSpace::FunctionSpace<MeshType, BasisFunctionType> FunctionSpaceType;
  
  std::array<dof_no_t, 9> dofIndices1;
  std::array<dof_no_t, 9> dofIndices1Reference = {0, 3, 6, 9, 12, 15, 18, 21, 24};
  FunctionSpaceType::getFaceDofs(Mesh::face_t::face0Minus, dofIndices1);  // 0-
  ASSERT_EQ(dofIndices1, dofIndices1Reference);
  
  std::array<dof_no_t, 9> dofIndices2;
  std::array<dof_no_t, 9> dofIndices2Reference = {2, 5, 8, 11, 14, 17, 20, 23, 26};
  FunctionSpaceType::getFaceDofs(Mesh::face_t::face0Plus, dofIndices2);  // 0+
  ASSERT_EQ(dofIndices2, dofIndices2Reference);
  
  std::array<dof_no_t, 9> dofIndices3;
  std::array<dof_no_t, 9> dofIndices3Reference = {0, 1, 2, 9, 10, 11, 18, 19, 20};
  FunctionSpaceType::getFaceDofs(Mesh::face_t::face1Minus, dofIndices3);  // 1-
  ASSERT_EQ(dofIndices3, dofIndices3Reference);
  
  std::array<dof_no_t, 9> dofIndices4;
  std::array<dof_no_t, 9> dofIndices4Reference = {6, 7, 8, 15, 16, 17, 24, 25, 26};
  FunctionSpaceType::getFaceDofs(Mesh::face_t::face1Plus, dofIndices4);  // 1+
  ASSERT_EQ(dofIndices4, dofIndices4Reference);
  
  std::array<dof_no_t, 9> dofIndices5;
  std::array<dof_no_t, 9> dofIndices5Reference = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  FunctionSpaceType::getFaceDofs(Mesh::face_t::face2Minus, dofIndices5);  // 2-
  ASSERT_EQ(dofIndices5, dofIndices5Reference);
  
  std::array<dof_no_t, 9> dofIndices6;
  std::array<dof_no_t, 9> dofIndices6Reference = {18, 19, 20, 21, 22, 23, 24, 25, 26};
  FunctionSpaceType::getFaceDofs(Mesh::face_t::face2Plus, dofIndices6);  // 2+
  ASSERT_EQ(dofIndices6, dofIndices6Reference);
  
}

TEST(FaceTest, normals3DCubeElement)
{
  typedef Mesh::StructuredDeformableOfDimension<3> MeshType;
  typedef BasisFunction::LagrangeOfOrder<1> BasisFunctionType;
  typedef FunctionSpace::FunctionSpace<MeshType, BasisFunctionType> FunctionSpaceType;
  
  
  std::string pythonConfig = R"(
# 3D
  
config = {
  "FiniteElementMethod" : {
    "nElements": 1,
    "nodePositions": [[0,0,0], [3,0,0], 
                      [0,3,0], [3,3,0], 
                      [0,0,3], [3,0,3],
                      [0,3,3], [3,3,3]],  # 2x2x2 nodes, 1 element
  }
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  FiniteElementMethod<
    MeshType,
    BasisFunctionType,
    Quadrature::Gauss<2>,
    Equation::None
  > problem(settings);
  
  std::shared_ptr<FunctionSpaceType> functionSpace = problem.functionSpace();
  
  {
    Vec3 xi({0.0, 0.5, 0.5});
    Vec3 normal1 = functionSpace->getNormal(Mesh::face_t::face0Minus, 0, xi);
    Vec3 normal1Reference({-1.0, 0.0, 0.0});
    for (int i = 0; i < 3; i++)
    {
      ASSERT_DOUBLE_EQ(normal1[i], normal1Reference[i]);
    }
    
    xi=Vec3({0.2, 0.4, 0.8});
    Vec3 normal2 = functionSpace->getNormal(Mesh::face_t::face0Minus, 0, xi);
    Vec3 normal2Reference({-1.0, 0.0, 0.0});
    for (int i = 0; i < 3; i++)
    {
      ASSERT_DOUBLE_EQ(normal2[i], normal2Reference[i]);
    }
  }
  
  {
    Vec3 xi({1.0, 0.5, 0.5});
    Vec3 normal1 = functionSpace->getNormal(Mesh::face_t::face0Plus, 0, xi);
    Vec3 normal1Reference({1.0, 0.0, 0.0});
    for (int i = 0; i < 3; i++)
    {
      ASSERT_DOUBLE_EQ(normal1[i], normal1Reference[i]);
    }
    
    xi=Vec3({0.2, 0.4, 0.8});
    Vec3 normal2 = functionSpace->getNormal(Mesh::face_t::face0Plus, 0, xi);
    Vec3 normal2Reference({1.0, 0.0, 0.0});
    for (int i = 0; i < 3; i++)
    {
      ASSERT_DOUBLE_EQ(normal2[i], normal2Reference[i]);
    }
  }
  
  {
    Vec3 xi({0.5, 0.0, 0.5});
    Vec3 normal1 = functionSpace->getNormal(Mesh::face_t::face1Minus, 0, xi);
    Vec3 normal1Reference({0.0, -1.0, 0.0});
    for (int i = 0; i < 3; i++)
    {
      ASSERT_DOUBLE_EQ(normal1[i], normal1Reference[i]);
    }
    
    xi=Vec3({0.2, 0.4, 0.8});
    Vec3 normal2 = functionSpace->getNormal(Mesh::face_t::face1Minus, 0, xi);
    Vec3 normal2Reference({0.0, -1.0, 0.0});
    for (int i = 0; i < 3; i++)
    {
      ASSERT_DOUBLE_EQ(normal2[i], normal2Reference[i]);
    }
  }
  
  {
    Vec3 xi({0.5, 1.0, 0.5});
    Vec3 normal1 = functionSpace->getNormal(Mesh::face_t::face1Plus, 0, xi);
    Vec3 normal1Reference({0.0, 1.0, 0.0});
    for (int i = 0; i < 3; i++)
    {
      ASSERT_DOUBLE_EQ(normal1[i], normal1Reference[i]);
    }
    
    xi=Vec3({0.2, 0.4, 0.8});
    Vec3 normal2 = functionSpace->getNormal(Mesh::face_t::face1Plus, 0, xi);
    Vec3 normal2Reference({0.0, 1.0, 0.0});
    for (int i = 0; i < 3; i++)
    {
      ASSERT_DOUBLE_EQ(normal2[i], normal2Reference[i]);
    }
  }
  {
    Vec3 xi({0.5, 0.5, 0.0});
    Vec3 normal1 = functionSpace->getNormal(Mesh::face_t::face2Minus, 0, xi);
    Vec3 normal1Reference({0.0, 0.0, -1.0});
    for (int i = 0; i < 3; i++)
    {
      ASSERT_DOUBLE_EQ(normal1[i], normal1Reference[i]);
    }
    
    xi=Vec3({0.2, 0.4, 0.8});
    Vec3 normal2 = functionSpace->getNormal(Mesh::face_t::face2Minus, 0, xi);
    Vec3 normal2Reference({0.0, 0.0, -1.0});
    for (int i = 0; i < 3; i++)
    {
      ASSERT_DOUBLE_EQ(normal2[i], normal2Reference[i]);
    }
  }
  
  {
    Vec3 xi({0.5, 0.5, 1.0});
    Vec3 normal1 = functionSpace->getNormal(Mesh::face_t::face2Plus, 0, xi);
    Vec3 normal1Reference({0.0, 0.0, 1.0});
    for (int i = 0; i < 3; i++)
    {
      ASSERT_DOUBLE_EQ(normal1[i], normal1Reference[i]);
    }
    
    xi=Vec3({0.2, 0.4, 0.8});
    Vec3 normal2 = functionSpace->getNormal(Mesh::face_t::face2Plus, 0, xi);
    Vec3 normal2Reference({0.0, 0.0, 1.0});
    for (int i = 0; i < 3; i++)
    {
      ASSERT_DOUBLE_EQ(normal2[i], normal2Reference[i]);
    }
  }
}

TEST(FaceTest, normals3DDistortedElement)
{
  typedef Mesh::StructuredDeformableOfDimension<3> MeshType;
  typedef BasisFunction::LagrangeOfOrder<1> BasisFunctionType;
  typedef FunctionSpace::FunctionSpace<MeshType, BasisFunctionType> FunctionSpaceType;
  
  
  std::string pythonConfig = R"(
# 3D
  
config = {
  "FiniteElementMethod" : {
    "nElements": 1,
    "nodePositions": [[0,0,0], [2,0,0], 
                      [0,2,0], [2,2,0], 
                      [1,0,2], [2,0,2],
                      [1,2,2], [2,2,2]],  # 2x2x2 nodes, 1 element
  }
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  FiniteElementMethod<
    MeshType,
    BasisFunctionType,
    Quadrature::Gauss<2>,
    Equation::None
  > problem(settings);
  problem.initialize();
  std::shared_ptr<FunctionSpaceType> functionSpace = std::static_pointer_cast<FunctionSpaceType>(problem.functionSpace());
  
  {
    Vec3 xi({0.0, 0.5, 0.5});
    Vec3 normal1 = functionSpace->getNormal(Mesh::face_t::face0Minus, 0, xi);
    double angle = atan2(1.,-2.);
    Vec3 normal1Reference({cos(angle), 0.0, sin(angle)});
    MathUtility::normalize<3>(normal1Reference);
    
    LOG(DEBUG) << " angle: " << angle << ", normal: " << normal1;
    
    for (int i = 0; i < 3; i++)
    {
      ASSERT_DOUBLE_EQ(normal1[i], normal1Reference[i]);
    }
    
    xi=Vec3({0.5, 0.0, 0.0});
    Vec3 normal2 = functionSpace->getNormal(Mesh::face_t::face0Minus, 0, xi);
    double angle2 = M_PI - angle;
    double angle3 = M_PI - angle2/2.;
    
    LOG(DEBUG) << " angle2: " << angle2 << ", angle3: " << angle3;
    
    Vec3 normal2Reference({-2, 0.0, 0.5});
    MathUtility::normalize<3>(normal2Reference);
    
    LOG(DEBUG) << " normal: " << normal2 << " (angle: " << atan2(normal2[2], normal2[0]) << ")";
    for (int i = 0; i < 3; i++)
    {
      ASSERT_DOUBLE_EQ(normal2[i], normal2Reference[i]);
    }
  }
}

TEST(FaceTest, normals3DDistortedElement2)
{
  typedef Mesh::StructuredDeformableOfDimension<3> MeshType;
  typedef BasisFunction::LagrangeOfOrder<1> BasisFunctionType;
  typedef FunctionSpace::FunctionSpace<MeshType, BasisFunctionType> FunctionSpaceType;
  
  
  std::string pythonConfig = R"(
# 3D
  
config = {
  "FiniteElementMethod" : {
    "nElements": 1,
    "nodePositions": [[0,0,0], [2,0,0],
                      [1,3,0], [3,2,0], 
                      [0,0,2], [2,0,2],
                      [1,3,3.5], [3,2,3]],  # 2x2x2 nodes, 1 element
  }
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  FiniteElementMethod<
    MeshType,
    BasisFunctionType,
    Quadrature::Gauss<2>,
    Equation::None
  > problem(settings);
  
  std::shared_ptr<FunctionSpaceType> functionSpace = std::static_pointer_cast<FunctionSpaceType>(problem.functionSpace());
  
  {
    Vec3 xi({0.0, 0.5, 0.5});
    Vec3 normal1 = functionSpace->getNormal(Mesh::face_t::face0Minus, 0, xi);
    Vec3 normal1Reference({-3.,1.,0.});
    MathUtility::normalize<3>(normal1Reference);
    LOG(DEBUG) << " 0- normal: " << normal1;
    for (int i = 0; i < 3; i++)
    {
      ASSERT_DOUBLE_EQ(normal1[i], normal1Reference[i]);
    }
    
    xi=Vec3({1.0, 0.5, 0.5});
    Vec3 normal2 = functionSpace->getNormal(Mesh::face_t::face0Plus, 0, xi);
    Vec3 normal2Reference({2., -1.0, 0.0});
    MathUtility::normalize<3>(normal2Reference);
    
    LOG(DEBUG) << " 0+ normal: " << normal2;
    for (int i = 0; i < 3; i++)
    {
      ASSERT_DOUBLE_EQ(normal2[i], normal2Reference[i]);
    }
  }
  
  {
    Vec3 xi({0.5, 0.0, 0.5});
    Vec3 normal1 = functionSpace->getNormal(Mesh::face_t::face1Minus, 0, xi);
    Vec3 normal1Reference({0.,-1.,0.});
    MathUtility::normalize<3>(normal1Reference);
    LOG(DEBUG) << " 1- normal: " << normal1;
    for (int i = 0; i < 3; i++)
    {
      ASSERT_DOUBLE_EQ(normal1[i], normal1Reference[i]);
    }
    
    xi=Vec3({0.5, 1.0, 0.5});
    Vec3 normal2 = functionSpace->getNormal(Mesh::face_t::face1Plus, 0, xi);
    Vec3 normal2Reference({1.,2.,0.});
    MathUtility::normalize<3>(normal2Reference);
    
    LOG(DEBUG) << " 1+ normal: " << normal2;
    for (int i = 0; i < 3; i++)
    {
      ASSERT_DOUBLE_EQ(normal2[i], normal2Reference[i]);
    }
  }
  
  {
    Vec3 xi=Vec3({0.5, 0.5, 1.0});
    Vec3 normal2 = functionSpace->getNormal(Mesh::face_t::face2Plus, 0, xi);
    Vec3 normal2Reference({0.,-1.,2.});
    MathUtility::normalize<3>(normal2Reference);
    
    LOG(DEBUG) << " 2+ normal: " << normal2;
    for (int i = 0; i < 3; i++)
    {
      ASSERT_DOUBLE_EQ(normal2[i], normal2Reference[i]);
    }
  }
}
};

