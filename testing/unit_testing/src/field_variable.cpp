#include <Python.h>  // this has to be the first included header

#include <iostream>
#include <cstdlib>
#include <fstream>

#include "gtest/gtest.h"
#include "opendihu.h"
#include "arg.h"
#include "stiffness_matrix_tester.h"
#include "node_positions_tester.h"

namespace Testing
{
  
TEST(FieldVariableTest, StructuredDeformable)
{
  // explicit mesh with node positions
  std::string pythonConfig = R"(
# Laplace 1D
config = {
  "Meshes" : {
    "testMesh": {
      "nElements": [2,2],
      "physicalExtent": [1.0,1.0],
    }
  },
  "FiniteElementMethod" : {
    "relativeTolerance": 1e-15,
    "meshName": "testMesh",
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<2>,
    BasisFunction::LagrangeOfOrder<2>,
    Quadrature::Gauss<2>,
    Equation::Static::Laplace
  > finiteElementMethod(settings);
  
  typedef Mesh::StructuredDeformableOfDimension<2> MeshType;
  typedef BasisFunction::LagrangeOfOrder<2> BasisFunctionType;
  
  typedef BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType> BasisOnMeshType;
  typedef std::shared_ptr<FieldVariable::FieldVariableBase<BasisOnMeshType>> FieldVariableBaseType;
  typedef std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,1>> FieldVariable1Type;
  typedef std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,2>> FieldVariable2Type;
  typedef std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,3>> FieldVariable3Type;
  
  std::shared_ptr<BasisOnMeshType> mesh = std::static_pointer_cast<BasisOnMeshType>(finiteElementMethod.mesh());
  
  // 5x5 nodes, 2 dofs/node, 18 dofs/element, 50 dofs
  
  FieldVariableBaseType aBase = mesh->createFieldVariable("a", {"x","y"});  
  FieldVariableBaseType bBase = mesh->createFieldVariable("b", 2);
  FieldVariableBaseType cBase = mesh->createFieldVariable("c");
  FieldVariable3Type d = mesh->template createFieldVariable<3>("d");
  
  
  FieldVariable2Type a = std::static_pointer_cast<FieldVariable::FieldVariable<BasisOnMeshType,2>>(aBase);
  FieldVariable2Type b = std::static_pointer_cast<FieldVariable::FieldVariable<BasisOnMeshType,2>>(bBase);
  FieldVariable1Type c = std::static_pointer_cast<FieldVariable::FieldVariable<BasisOnMeshType,1>>(cBase);
  
  // set all to 0.0
  a->setValues(0.0);
  a->finishVectorManipulation();
  
  ASSERT_EQ(a->getValue(0, 1), 0.0);
  ASSERT_EQ(a->getValue(0, 5), 0.0);
  ASSERT_EQ(a->getValue(1, 8), 0.0);
  
  // set all components, 1 dof
  Vec2 v0 = {1.0, 2.0};
  a->setValue(5, v0);
  a->finishVectorManipulation();
  
  ASSERT_EQ(a->getValue(0, 5), v0[0]);
  ASSERT_EQ(a->getValue(1, 5), v0[1]);
  
  // set all components, multiple dofs
  std::vector<Vec2> v1 = {Vec2({1.0, 2.0}), Vec2({3.0, 4.0}), Vec2({5.0, 6.0}), Vec2({7.0, 8.0}), Vec2({9.0, 10.0})};
  std::vector<dof_no_t> dofs = {4,8,2,7,10};
  a->setValues(dofs, v1);
  a->finishVectorManipulation();
  
  ASSERT_EQ(a->getValue(0,4), 1.0);
  ASSERT_EQ(a->getValue(1,4), 2.0);
  ASSERT_EQ(a->getValue(0,8), 3.0);
  ASSERT_EQ(a->getValue(1,8), 4.0);
  ASSERT_EQ(a->getValue(0,2), 5.0);
  ASSERT_EQ(a->getValue(1,2), 6.0);
  ASSERT_EQ(a->getValue(0,7), 7.0);
  ASSERT_EQ(a->getValue(1,7), 8.0);
  ASSERT_EQ(a->getValue(0,10), 9.0);
  ASSERT_EQ(a->getValue(1,10), 10.0);

  // set all to 0.0
  c->setValues(0.0);
  c->finishVectorManipulation();
  
  // set single component, 1 dof
  c->setValue(5, 5.0);
  c->finishVectorManipulation();
  ASSERT_EQ(c->getValue(0,5), 5.0);
  ASSERT_EQ(c->getValue(5), 5.0);
  ASSERT_EQ(c->getValue(0), 0.0);
  
  // set single component, multiple dofs
  std::vector<double> v2 = {1.0, 2.0, 3.0, 4.0, 5.0};
  c->setValues(dofs, v2);
  c->finishVectorManipulation();
  
  ASSERT_EQ(c->getValue(4), 1.0);
  ASSERT_EQ(c->getValue(8), 2.0);
  ASSERT_EQ(c->getValue(2), 3.0);
  ASSERT_EQ(c->getValue(7), 4.0);
  ASSERT_EQ(c->getValue(10), 5.0);
  
  std::vector<double> values10 = {11., 12., 13., 14.};
  std::vector<dof_no_t> dofs10 = {11, 12, 13, 14};
  c->setValues(dofs10, values10);
  c->finishVectorManipulation();
  
  /* values of c
     0.0, 0.0, 3.0, 0.0, 1.0,
     5.0, 0.0, 4.0, 2.0, 0.0,
     5.0,11.0,12.0,13.0,14.0,
     0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,
   */
// dofs element 1: 2,3,4, 7,8,9, 12,13,14
  
  std::array<double,BasisOnMeshType::nDofsPerElement()> values8;
  c->getElementValues(1.0, values8);
  
  std::array<double,BasisOnMeshType::nDofsPerElement()> reference8 = {
    3.0, 0.0, 1.0,  4.0, 2.0, 0.0,  12.0, 13.0, 14.0
  };
  ASSERT_EQ(values8, reference8);
  
  
  /* values of a:
    0.0, 0.0,  0.0, 0.0,  5.0, 6.0,  0.0, 0.0,  1.0, 2.0,
    1.0, 2.0,  0.0, 0.0,  7.0, 8.0,  3.0, 4.0,  0.0, 0.0,
    9.0,10.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,
    0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,
    0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0
    */
  
  // for a component get all values
  std::vector<double> values;
  a->getValues(0, values);
  
  std::vector<double> reference0 = {
    0.0, 0.0, 5.0, 0.0, 1.0,
    1.0, 0.0, 7.0, 3.0, 0.0,
    9.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0
  };
  ASSERT_EQ(values, reference0);
  
  values.clear();
  a->getValues(1, values);
  
  std::vector<double> reference1 = {
     0.0, 0.0, 6.0, 0.0, 2.0,
     2.0, 0.0, 8.0, 4.0, 0.0,
    10.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0
  };
  ASSERT_EQ(values, reference1);
  
  
  // one component, multiple dofs 
  std::array<dof_no_t,3> dofs2 = {1,7,2};
  std::array<double,3> values2;
  
  a->template getValues<3>(0,dofs2,values2);
  std::array<double,3> reference2{0.0,7.0,5.0};
  ASSERT_EQ(values2,reference2);
  
  a->template getValues<3>(1,dofs2,values2);
  std::array<double,3> reference3{0.0,8.0,6.0};
  ASSERT_EQ(values2,reference3);
  
  // same but with vectors
  std::vector<dof_no_t> dofs3 = {1,7,2};
  std::vector<double> values3 = {-1.0,-1.0};
  
  a->getValues(0,dofs3,values3);
  std::vector<double> rreference3 = {-1.0,-1.0,0.0,7.0,5.0};
  ASSERT_EQ(values3,rreference3);
  
  a->getValues(1,dofs3,values3);
  std::vector<double> rreference33 = {-1.0,-1.0,0.0,7.0,5.0,0.0,8.0,6.0};
  ASSERT_EQ(values3,rreference33);
  
  // all components, multiple dofs
  std::array<Vec2,3> values4;
  a->template getValues<3>(dofs2, values4);
  std::array<Vec2,3> reference4 = {Vec2({0.0,0.0}),Vec2({7.0,8.0}),Vec2({5.0,6.0})};
  ASSERT_EQ(values4,reference4);
  
  
//  node global nos.
//  10__11_12__13_14
//  ++==+==++==+==++
// 5||_6|_7||_8|__||9
//  ||__|__||__|__||
//   0  1  2   3  4
// dofs element 1: 2,3,4, 7,8,9, 12,13,14
  
  //! for a specific component, get the values corresponding to all element-local dofs
  // component 0
  std::array<double,BasisOnMeshType::nDofsPerElement()> values5;
  a->getElementValues(0, 1, values5);
  std::array<double,BasisOnMeshType::nDofsPerElement()> reference5 = {5.0, 0.0, 1.0, 7.0, 3.0, 0.0, 0.0, 0.0, 0.0};
  ASSERT_EQ(values5,reference5);
  
  // component 1
  std::array<double,BasisOnMeshType::nDofsPerElement()> values6;
  a->getElementValues(1, 1, values6);
  std::array<double,BasisOnMeshType::nDofsPerElement()> reference6 = {6.0, 0.0, 2.0, 8.0, 4.0, 0.0, 0.0, 0.0, 0.0};
  ASSERT_EQ(values6,reference6);
  
  
  
  //! get the values corresponding to all element-local dofs for all components
  std::array<Vec2,BasisOnMeshType::nDofsPerElement()> values7;
  a->getElementValues(1, values7);
  std::array<Vec2,BasisOnMeshType::nDofsPerElement()> reference7 = {Vec2({5.0, 6.0}), Vec2({0.0, 0.0}), Vec2({1.0, 2.0}), Vec2({7.0, 8.0}), Vec2({3.0, 4.0}), Vec2({0.0, 0.0}), Vec2({0.0, 0.0}), Vec2({0.0, 0.0}), Vec2({0.0, 0.0})};
  
  ASSERT_EQ(values7,reference7);
  
  // copy a to b
  b->setValues(*a);
  b->finishVectorManipulation();
  b->getElementValues(1, values7);
  ASSERT_EQ(values7,reference7);
  
  /*
   * 
  //! copy the values from another field variable of the same type
  void setValues(FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents> &rhs);
  */
  
  /*
  //! create a non-geometry field field variable with no values being set, with given component names
  std::shared_ptr<FieldVariable::FieldVariableBase<BasisOnMesh<MeshType,BasisFunctionType>>> createFieldVariable(std::string name, std::vector<std::string> componentNames);
  
  //! create a non-geometry field field variable with no values being set, with given number of components, the component names will be the numbers
  std::shared_ptr<FieldVariable::FieldVariableBase<BasisOnMesh<MeshType,BasisFunctionType>>> createFieldVariable(std::string name, int nComponents=1);
  
  //! create a non-geometry field field variable with no values being set, with given number of components, the component names will be the numbers
  template <int nComponents>
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMesh<MeshType,BasisFunctionType>,nComponents>> createFieldVariable(std::string name);
  */
  
  /*
   * 
  //! set values for all components for dofs, after all calls to setValue(s), finishVectorManipulation has to be called to apply the cached changes
  void setValues(std::vector<dof_no_t> &dofGlobalNos, std::vector<std::array<double,nComponents>> &values, InsertMode petscInsertMode=INSERT_VALUES)
  
  //! set a single dof (all components) , after all calls to setValue(s), finishVectorManipulation has to be called to apply the cached changes
  void setValue(dof_no_t dofGlobalNo, std::array<double,nComponents> &value, InsertMode petscInsertMode=INSERT_VALUES)
 
   * 
  //! for a specific component, get all values
  //! @param onlyNodalValues: if this is true, for Hermite only the non-derivative values are retrieved
  void getValues(int componentNo, std::vector<double> &values, bool onlyNodalValues=false);
  
  //! for a specific component, get values from their global dof no.s, as array, therefore templated by the number of elements, N, to retrieve
  template<int N>
  void getValues(int componentNo, std::array<dof_no_t,N> dofGlobalNo, std::array<double,N> &values);
  
  //! for a specific component, get values from their global dof no.s, as vector
  void getValues(int componentNo, std::vector<dof_no_t> dofGlobalNo, std::vector<double> &values);
  
  //! get values from their global dof no.s for all components
  template<int N>
  void getValues(std::array<dof_no_t,N> dofGlobalNo, std::array<std::array<double,nComponents>,N> &values)
  
  //! for a specific component, get the values corresponding to all element-local dofs
  template<int N>
  void getElementValues(int componentNo, element_no_t elementNo, std::array<double,BasisOnMeshType::nDofsPerElement()> &values);
  
  //! get the values corresponding to all element-local dofs for all components
  void getElementValues(element_no_t elementNo, std::array<std::array<double,nComponents>,BasisOnMeshType::nDofsPerElement()> &values)
 
  //! for a specific component, get a single value from global dof no.
  double getValue(int componentNo, node_no_t dofGlobalNo);




  std::array<double,nComponents> getValue(node_no_t dofGlobalNo);
  
  
  
  
  
  //! get the values corresponding to all element-local dofs for all components
  void getElementValues(element_no_t elementNo, std::array<double,BasisOnMeshType::nDofsPerElement()> &values)
  
  //! get a single value from global dof no. for all components
  double getValue(node_no_t dofGlobalNo);
  
  //! set a single dof (all components) , after all calls to setValue(s), finishVectorManipulation has to be called to apply the cached changes
  void setValue(dof_no_t dofGlobalNo, double value, InsertMode petscInsertMode=INSERT_VALUES)

  //! set values for all components for dofs, after all calls to setValue(s), finishVectorManipulation has to be called to apply the cached changes
  void setValues(std::vector<dof_no_t> &dofGlobalNos, std::vector<double> &values, InsertMode petscInsertMode=INSERT_VALUES)
*/
}

TEST(FieldVariableTest, StructuredRegularFixed)
{
  // explicit mesh with node positions
  std::string pythonConfig = R"(
# Laplace 1D
config = {
  "Meshes" : {
    "testMesh": {
      "nElements": [2,2],
      "physicalExtent": [1.0,1.0],
    }
  },
  "FiniteElementMethod" : {
    "relativeTolerance": 1e-15,
    "meshName": "testMesh",
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredRegularFixedOfDimension<2>,
    BasisFunction::LagrangeOfOrder<2>,
    Quadrature::Gauss<2>,
    Equation::Static::Laplace
  > finiteElementMethod(settings);
  
  typedef Mesh::StructuredRegularFixedOfDimension<2> MeshType;
  typedef BasisFunction::LagrangeOfOrder<2> BasisFunctionType;
  
  typedef BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType> BasisOnMeshType;
  typedef std::shared_ptr<FieldVariable::FieldVariableBase<BasisOnMeshType>> FieldVariableBaseType;
  typedef std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,1>> FieldVariable1Type;
  typedef std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,2>> FieldVariable2Type;
  typedef std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,3>> FieldVariable3Type;
  
  std::shared_ptr<BasisOnMeshType> mesh = std::static_pointer_cast<BasisOnMeshType>(finiteElementMethod.mesh());
  
  // 5x5 nodes, 2 dofs/node, 18 dofs/element, 50 dofs
  
  FieldVariableBaseType aBase = mesh->createFieldVariable("a", {"x","y"});  
  FieldVariableBaseType bBase = mesh->createFieldVariable("b", 2);
  FieldVariableBaseType cBase = mesh->createFieldVariable("c");
  FieldVariable3Type d = mesh->template createFieldVariable<3>("d");
  
  
  FieldVariable2Type a = std::static_pointer_cast<FieldVariable::FieldVariable<BasisOnMeshType,2>>(aBase);
  FieldVariable2Type b = std::static_pointer_cast<FieldVariable::FieldVariable<BasisOnMeshType,2>>(bBase);
  FieldVariable1Type c = std::static_pointer_cast<FieldVariable::FieldVariable<BasisOnMeshType,1>>(cBase);
  
  // set all to 0.0
  a->setValues(0.0);
  a->finishVectorManipulation();
  
  ASSERT_EQ(a->getValue(0, 1), 0.0);
  ASSERT_EQ(a->getValue(0, 5), 0.0);
  ASSERT_EQ(a->getValue(1, 8), 0.0);
  
  // set all components, 1 dof
  Vec2 v0 = {1.0, 2.0};
  a->setValue(5, v0);
  a->finishVectorManipulation();
  
  ASSERT_EQ(a->getValue(0, 5), v0[0]);
  ASSERT_EQ(a->getValue(1, 5), v0[1]);
  
  // set all components, multiple dofs
  std::vector<Vec2> v1 = {Vec2({1.0, 2.0}), Vec2({3.0, 4.0}), Vec2({5.0, 6.0}), Vec2({7.0, 8.0}), Vec2({9.0, 10.0})};
  std::vector<dof_no_t> dofs = {4,8,2,7,10};
  a->setValues(dofs, v1);
  a->finishVectorManipulation();
  
  ASSERT_EQ(a->getValue(0,4), 1.0);
  ASSERT_EQ(a->getValue(1,4), 2.0);
  ASSERT_EQ(a->getValue(0,8), 3.0);
  ASSERT_EQ(a->getValue(1,8), 4.0);
  ASSERT_EQ(a->getValue(0,2), 5.0);
  ASSERT_EQ(a->getValue(1,2), 6.0);
  ASSERT_EQ(a->getValue(0,7), 7.0);
  ASSERT_EQ(a->getValue(1,7), 8.0);
  ASSERT_EQ(a->getValue(0,10), 9.0);
  ASSERT_EQ(a->getValue(1,10), 10.0);

  // set all to 0.0
  c->setValues(0.0);
  c->finishVectorManipulation();
  
  // set single component, 1 dof
  c->setValue(5, 5.0);
  c->finishVectorManipulation();
  ASSERT_EQ(c->getValue(0,5), 5.0);
  ASSERT_EQ(c->getValue(5), 5.0);
  ASSERT_EQ(c->getValue(0), 0.0);
  
  // set single component, multiple dofs
  std::vector<double> v2 = {1.0, 2.0, 3.0, 4.0, 5.0};
  c->setValues(dofs, v2);
  c->finishVectorManipulation();
  
  ASSERT_EQ(c->getValue(4), 1.0);
  ASSERT_EQ(c->getValue(8), 2.0);
  ASSERT_EQ(c->getValue(2), 3.0);
  ASSERT_EQ(c->getValue(7), 4.0);
  ASSERT_EQ(c->getValue(10), 5.0);
  
  std::vector<double> values10 = {11., 12., 13., 14.};
  std::vector<dof_no_t> dofs10 = {11, 12, 13, 14};
  c->setValues(dofs10, values10);
  c->finishVectorManipulation();
  
  /* values of c
     0.0, 0.0, 3.0, 0.0, 1.0,
     5.0, 0.0, 4.0, 2.0, 0.0,
     5.0,11.0,12.0,13.0,14.0,
     0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,
   */
// dofs element 1: 2,3,4, 7,8,9, 12,13,14
  
  std::array<double,BasisOnMeshType::nDofsPerElement()> values8;
  c->getElementValues(1.0, values8);
  
  std::array<double,BasisOnMeshType::nDofsPerElement()> reference8 = {
    3.0, 0.0, 1.0,  4.0, 2.0, 0.0,  12.0, 13.0, 14.0
  };
  ASSERT_EQ(values8, reference8);
  
  
  /* values of a:
    0.0, 0.0,  0.0, 0.0,  5.0, 6.0,  0.0, 0.0,  1.0, 2.0,
    1.0, 2.0,  0.0, 0.0,  7.0, 8.0,  3.0, 4.0,  0.0, 0.0,
    9.0,10.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,
    0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,
    0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0
    */
  
  // for a component get all values
  std::vector<double> values;
  a->getValues(0, values);
  
  std::vector<double> reference0 = {
    0.0, 0.0, 5.0, 0.0, 1.0,
    1.0, 0.0, 7.0, 3.0, 0.0,
    9.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0
  };
  ASSERT_EQ(values, reference0);
  
  values.clear();
  a->getValues(1, values);
  
  std::vector<double> reference1 = {
     0.0, 0.0, 6.0, 0.0, 2.0,
     2.0, 0.0, 8.0, 4.0, 0.0,
    10.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0
  };
  ASSERT_EQ(values, reference1);
  
  
  // one component, multiple dofs 
  std::array<dof_no_t,3> dofs2 = {1,7,2};
  std::array<double,3> values2;
  
  a->template getValues<3>(0,dofs2,values2);
  std::array<double,3> reference2{0.0,7.0,5.0};
  ASSERT_EQ(values2,reference2);
  
  a->template getValues<3>(1,dofs2,values2);
  std::array<double,3> reference3{0.0,8.0,6.0};
  ASSERT_EQ(values2,reference3);
  
  // same but with vectors
  std::vector<dof_no_t> dofs3 = {1,7,2};
  std::vector<double> values3 = {-1.0,-1.0};
  
  a->getValues(0,dofs3,values3);
  std::vector<double> rreference3 = {-1.0,-1.0,0.0,7.0,5.0};
  ASSERT_EQ(values3,rreference3);
  
  a->getValues(1,dofs3,values3);
  std::vector<double> rreference33 = {-1.0,-1.0,0.0,7.0,5.0,0.0,8.0,6.0};
  ASSERT_EQ(values3,rreference33);
  
  // all components, multiple dofs
  std::array<Vec2,3> values4;
  a->template getValues<3>(dofs2, values4);
  std::array<Vec2,3> reference4 = {Vec2({0.0,0.0}),Vec2({7.0,8.0}),Vec2({5.0,6.0})};
  ASSERT_EQ(values4,reference4);
  
  
//  node global nos.
//  10__11_12__13_14
//  ++==+==++==+==++
// 5||_6|_7||_8|__||9
//  ||__|__||__|__||
//   0  1  2   3  4
// dofs element 1: 2,3,4, 7,8,9, 12,13,14
  
  //! for a specific component, get the values corresponding to all element-local dofs
  // component 0
  std::array<double,BasisOnMeshType::nDofsPerElement()> values5;
  a->getElementValues(0, 1, values5);
  std::array<double,BasisOnMeshType::nDofsPerElement()> reference5 = {5.0, 0.0, 1.0, 7.0, 3.0, 0.0, 0.0, 0.0, 0.0};
  ASSERT_EQ(values5,reference5);
  
  // component 1
  std::array<double,BasisOnMeshType::nDofsPerElement()> values6;
  a->getElementValues(1, 1, values6);
  std::array<double,BasisOnMeshType::nDofsPerElement()> reference6 = {6.0, 0.0, 2.0, 8.0, 4.0, 0.0, 0.0, 0.0, 0.0};
  ASSERT_EQ(values6,reference6);
  
  
  
  //! get the values corresponding to all element-local dofs for all components
  std::array<Vec2,BasisOnMeshType::nDofsPerElement()> values7;
  a->getElementValues(1, values7);
  std::array<Vec2,BasisOnMeshType::nDofsPerElement()> reference7 = {Vec2({5.0, 6.0}), Vec2({0.0, 0.0}), Vec2({1.0, 2.0}), Vec2({7.0, 8.0}), Vec2({3.0, 4.0}), Vec2({0.0, 0.0}), Vec2({0.0, 0.0}), Vec2({0.0, 0.0}), Vec2({0.0, 0.0})};
  
  ASSERT_EQ(values7,reference7);
  
  // copy a to b
  b->setValues(*a);
  b->finishVectorManipulation();
  b->getElementValues(1, values7);
  ASSERT_EQ(values7,reference7);
  
}

TEST(FieldVariableTest, UnstructuredDeformable)
{
  // explicit mesh with node positions
  std::string pythonConfig = R"(
# Laplace 1D
import numpy as np

nodePositions = []
for y in np.linspace(0,1,5):
  for x in np.linspace(0,1,5):
    nodePositions.append([x,y])
    
# node numbers
# 0 1 2 3 4   
# 5 6 7 8 9    
# 0 1 2 3 4
# 5 6 7 8 9
# 0 1 2 3 4
    
config = {
  "Meshes" : {
    "testMesh": {
      "nElements": 4,
      "nodePositions": nodePositions,
      "elements": [[0,1,2,5,6,7,10,11,12], [2,3,4,7,8,9,12,13,14], [10,11,12,15,16,17,20,21,22], [12,13,14,17,18,19,22,23,24]],
    }
  },
  "FiniteElementMethod" : {
    "relativeTolerance": 1e-15,
    "meshName": "testMesh",
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);
  
  SpatialDiscretization::FiniteElementMethod<
    Mesh::UnstructuredDeformableOfDimension<2>,
    BasisFunction::LagrangeOfOrder<2>,
    Quadrature::Gauss<2>,
    Equation::Static::Laplace
  > finiteElementMethod(settings);
  
  typedef Mesh::UnstructuredDeformableOfDimension<2> MeshType;
  typedef BasisFunction::LagrangeOfOrder<2> BasisFunctionType;
  
  typedef BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType> BasisOnMeshType;
  typedef std::shared_ptr<FieldVariable::FieldVariableBase<BasisOnMeshType>> FieldVariableBaseType;
  typedef std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,1>> FieldVariable1Type;
  typedef std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,2>> FieldVariable2Type;
  typedef std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,3>> FieldVariable3Type;
  
  std::shared_ptr<BasisOnMeshType> mesh = std::static_pointer_cast<BasisOnMeshType>(finiteElementMethod.mesh());
  
  // 5x5 nodes, 2 dofs/node, 18 dofs/element, 50 dofs
  
  FieldVariableBaseType aBase = mesh->createFieldVariable("a", {"x","y"});  
  FieldVariableBaseType bBase = mesh->createFieldVariable("b", 2);
  FieldVariableBaseType cBase = mesh->createFieldVariable("c");
  FieldVariable3Type d = mesh->template createFieldVariable<3>("d");
  
  
  FieldVariable2Type a = std::static_pointer_cast<FieldVariable::FieldVariable<BasisOnMeshType,2>>(aBase);
  FieldVariable2Type b = std::static_pointer_cast<FieldVariable::FieldVariable<BasisOnMeshType,2>>(bBase);
  FieldVariable1Type c = std::static_pointer_cast<FieldVariable::FieldVariable<BasisOnMeshType,1>>(cBase);
  
  // set all to 0.0
  a->setValues(0.0);
  a->finishVectorManipulation();
  
  ASSERT_EQ(a->getValue(0, 1), 0.0);
  ASSERT_EQ(a->getValue(0, 5), 0.0);
  ASSERT_EQ(a->getValue(1, 8), 0.0);
  
  // set all components, 1 dof
  Vec2 v0 = {1.0, 2.0};
  a->setValue(5, v0);
  a->finishVectorManipulation();
  
  ASSERT_EQ(a->getValue(0, 5), v0[0]);
  ASSERT_EQ(a->getValue(1, 5), v0[1]);
  
  // set all components, multiple dofs
  std::vector<Vec2> v1 = {Vec2({1.0, 2.0}), Vec2({3.0, 4.0}), Vec2({5.0, 6.0}), Vec2({7.0, 8.0}), Vec2({9.0, 10.0})};
  std::vector<dof_no_t> dofs = {4,8,2,7,10};
  a->setValues(dofs, v1);
  a->finishVectorManipulation();
  
  ASSERT_EQ(a->getValue(0,4), 1.0);
  ASSERT_EQ(a->getValue(1,4), 2.0);
  ASSERT_EQ(a->getValue(0,8), 3.0);
  ASSERT_EQ(a->getValue(1,8), 4.0);
  ASSERT_EQ(a->getValue(0,2), 5.0);
  ASSERT_EQ(a->getValue(1,2), 6.0);
  ASSERT_EQ(a->getValue(0,7), 7.0);
  ASSERT_EQ(a->getValue(1,7), 8.0);
  ASSERT_EQ(a->getValue(0,10), 9.0);
  ASSERT_EQ(a->getValue(1,10), 10.0);

  // set all to 0.0
  c->setValues(0.0);
  c->finishVectorManipulation();
  
  // set single component, 1 dof
  c->setValue(5, 5.0);
  c->finishVectorManipulation();
  ASSERT_EQ(c->getValue(0,5), 5.0);
  ASSERT_EQ(c->getValue(5), 5.0);
  ASSERT_EQ(c->getValue(0), 0.0);
  
  // set single component, multiple dofs
  std::vector<double> v2 = {1.0, 2.0, 3.0, 4.0, 5.0};
  c->setValues(dofs, v2);
  c->finishVectorManipulation();
  
  ASSERT_EQ(c->getValue(4), 1.0);
  ASSERT_EQ(c->getValue(8), 2.0);
  ASSERT_EQ(c->getValue(2), 3.0);
  ASSERT_EQ(c->getValue(7), 4.0);
  ASSERT_EQ(c->getValue(10), 5.0);
  
  std::vector<double> values10 = {11., 12., 13., 14.};
  std::vector<dof_no_t> dofs10 = {11, 12, 13, 14};
  c->setValues(dofs10, values10);
  c->finishVectorManipulation();
  
  /* values of c
     0.0, 0.0, 3.0, 0.0, 1.0,
     5.0, 0.0, 4.0, 2.0, 0.0,
     5.0,11.0,12.0,13.0,14.0,
     0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,
   */
// dofs element 1: 2,9,10, 5,11,12, 8,13,14
  
  std::array<double,BasisOnMeshType::nDofsPerElement()> values8;
  c->getElementValues(1.0, values8);
  
  std::array<double,BasisOnMeshType::nDofsPerElement()> reference8 = {
    3.0, 0.0, 5.0,  5.0,11.0,12.0,  2.0, 13.0, 14.0
  };
  ASSERT_EQ(values8, reference8);
  
  
  /* values of a:
    0.0, 0.0,  0.0, 0.0,  5.0, 6.0,  0.0, 0.0,  1.0, 2.0,
    1.0, 2.0,  0.0, 0.0,  7.0, 8.0,  3.0, 4.0,  0.0, 0.0,
    9.0,10.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,
    0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,
    0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0,  0.0, 0.0
    */
  
  // for a component get all values
  std::vector<double> values;
  a->getValues(0, values);
  
  std::vector<double> reference0 = {
    0.0, 0.0, 5.0, 0.0, 1.0,
    1.0, 0.0, 7.0, 3.0, 0.0,
    9.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0
  };
  ASSERT_EQ(values, reference0);
  
  values.clear();
  a->getValues(1, values);
  
  std::vector<double> reference1 = {
     0.0, 0.0, 6.0, 0.0, 2.0,
     2.0, 0.0, 8.0, 4.0, 0.0,
    10.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0
  };
  ASSERT_EQ(values, reference1);
  
  
  // one component, multiple dofs 
  std::array<dof_no_t,3> dofs2 = {1,7,2};
  std::array<double,3> values2;
  
  a->template getValues<3>(0,dofs2,values2);
  std::array<double,3> reference2{0.0,7.0,5.0};
  ASSERT_EQ(values2,reference2);
  
  a->template getValues<3>(1,dofs2,values2);
  std::array<double,3> reference3{0.0,8.0,6.0};
  ASSERT_EQ(values2,reference3);
  
  // same but with vectors
  std::vector<dof_no_t> dofs3 = {1,7,2};
  std::vector<double> values3 = {-1.0,-1.0};
  
  a->getValues(0,dofs3,values3);
  std::vector<double> rreference3 = {-1.0,-1.0,0.0,7.0,5.0};
  ASSERT_EQ(values3,rreference3);
  
  a->getValues(1,dofs3,values3);
  std::vector<double> rreference33 = {-1.0,-1.0,0.0,7.0,5.0,0.0,8.0,6.0};
  ASSERT_EQ(values3,rreference33);
  
  // all components, multiple dofs
  std::array<Vec2,3> values4;
  a->template getValues<3>(dofs2, values4);
  std::array<Vec2,3> reference4 = {Vec2({0.0,0.0}),Vec2({7.0,8.0}),Vec2({5.0,6.0})};
  ASSERT_EQ(values4,reference4);
  
  
// dof nos.
  
// 18 19 20 23 24   
// 15 16 17 21 22    
//  6  7  8 13 14
//  3  4  5 11 12
//  0  1  2  9 10
// dofs element 1: 2,9,10, 5,11,12, 8,13,14
  
  //! for a specific component, get the values corresponding to all element-local dofs
  // component 0
  std::array<double,BasisOnMeshType::nDofsPerElement()> values5;
  a->getElementValues(0, 1, values5);
  std::array<double,BasisOnMeshType::nDofsPerElement()> reference5 = {5.0, 0.0, 9.0, 1.0, 0.0, 0.0, 3.0, 0.0, 0.0};
  ASSERT_EQ(values5,reference5);
  
  // component 1
  std::array<double,BasisOnMeshType::nDofsPerElement()> values6;
  a->getElementValues(1, 1, values6);
  std::array<double,BasisOnMeshType::nDofsPerElement()> reference6 = {6.0, 0.0,10.0, 2.0, 0.0, 0.0, 4.0, 0.0, 0.0};
  ASSERT_EQ(values6,reference6);
  
  
  
  //! get the values corresponding to all element-local dofs for all components
  std::array<Vec2,BasisOnMeshType::nDofsPerElement()> values7;
  a->getElementValues(1, values7);
  std::array<Vec2,BasisOnMeshType::nDofsPerElement()> reference7 
    = {Vec2({5.0, 6.0}), Vec2({0.0, 0.0}), Vec2({9.0, 10.0}), Vec2({1.0, 2.0}), Vec2({0.0, 0.0}), Vec2({0.0, 0.0}), 
       Vec2({3.0, 4.0}), Vec2({0.0, 0.0}), Vec2({0.0, 0.0})};
  
  ASSERT_EQ(values7,reference7);
  
  // copy a to b
  b->setValues(*a);
  b->finishVectorManipulation();
  b->getElementValues(1, values7);
  ASSERT_EQ(values7,reference7);
  
}

}; // namespace
