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

TEST(FieldVariableTest, GenericFieldVariable)
{
  std::string pythonConfig = R"(
config = {
  "FiniteElementMethod" : {
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<2>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::Gauss<2>,
    Equation::None
  > finiteElementMethod(settings);
  finiteElementMethod.run();


  // std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::Generic>> createGenericFieldVariable(std::shared_ptr<Partition::Manager> partitionManager, int nEntries, std::string name);

  // get functionSpace manager object that are stored in the DihuContext object
  std::shared_ptr<Mesh::Manager> meshManager = settings.meshManager();

  // create the field variable with name "test"
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::Generic,1>> fieldVariable = meshManager->createGenericFieldVariable(5, "test");

  // set all values to 0
  fieldVariable->zeroEntries();

  // set test[0] = 10, test[4] = 20
  std::vector<dof_no_t> dofNos{0, 4};
  std::vector<double> values{10.0, 20.0};
  fieldVariable->setValues(dofNos, values);

  // set test[2] = 4
  fieldVariable->setValue(2, 4);

  // add +10 to test[4]
  fieldVariable->setValue(4, 10, ADD_VALUES);

  // add +1 to all values by Petsc function
  PetscErrorCode ierr;
  ierr = VecShift(fieldVariable->valuesGlobal(), 1); CHKERRV(ierr);

  // print field variable
  LOG(DEBUG) << *fieldVariable;

  // get all values
  std::vector<double> vectorValues;
  fieldVariable->getValuesWithoutGhosts(vectorValues);

  // because this test case is run serially, we get all 5 values
  ASSERT_EQ(vectorValues[0], 11.0);
  ASSERT_EQ(vectorValues[1], 1.0);
  ASSERT_EQ(vectorValues[2], 5.0);
  ASSERT_EQ(vectorValues[3], 1.0);
  ASSERT_EQ(vectorValues[4], 31.0);

}

TEST(FieldVariableTest, NonSquareDenseMatrix)
{
  std::string pythonConfig = R"(
config = {
  "FiniteElementMethod" : {
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredRegularFixedOfDimension<2>,
    BasisFunction::LagrangeOfOrder<>,
    Quadrature::Gauss<2>,
    Equation::None
  > finiteElementMethod(settings);
  finiteElementMethod.run();

  typedef FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<2>,BasisFunction::LagrangeOfOrder<1>> FunctionSpace1;
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<2>,BasisFunction::LagrangeOfOrder<2>> FunctionSpace2;

  // get functionSpace manager object that are stored in the DihuContext object
  std::shared_ptr<Mesh::Manager> meshManager = settings.meshManager();

  std::array<int,2> nElements1({2,2});
  std::array<int,2> nElements2({3,3});
  std::array<double,2> physicalExtent({0});
  std::shared_ptr<FunctionSpace1> functionSpace1 = meshManager->createFunctionSpace<FunctionSpace1>("functionSpace1", nElements1, physicalExtent);
  std::shared_ptr<FunctionSpace2> functionSpace2 = meshManager->createFunctionSpace<FunctionSpace2>("functionSpace2", nElements2, physicalExtent);

  // create dense matrix, 9x49
  std::shared_ptr<PartitionedPetscMat<FunctionSpace1,FunctionSpace2>> matrix
    = std::make_shared<PartitionedPetscMat<FunctionSpace1,FunctionSpace2>>(functionSpace1->meshPartition(), functionSpace2->meshPartition(), 1, "T");

  // zero out all entries
  matrix->zeroEntries();

  // set some values in the matrix
  std::vector<PetscInt> rowIndices    = {0, 2};
  std::vector<PetscInt> columnIndices = {0, 48};
  std::vector<double> values(rowIndices.size() * columnIndices.size(), 1.0);
  matrix->setValues(rowIndices.size(), rowIndices.data(), columnIndices.size(), columnIndices.data(), values.data(), INSERT_VALUES);

  matrix->assembly(MAT_FINAL_ASSEMBLY);

  // get all values in the matrix
  std::vector<PetscInt> resultRowIndices(9);
  std::iota(resultRowIndices.begin(), resultRowIndices.end(), 0);
  std::vector<PetscInt> resultColumnIndices(49);
  std::iota(resultColumnIndices.begin(), resultColumnIndices.end(), 0);
  std::vector<double> resultValues(resultRowIndices.size() * resultColumnIndices.size());

  matrix->getValues(resultRowIndices.size(), resultRowIndices.data(), resultColumnIndices.size(), resultColumnIndices.data(), resultValues.data());

  LOG(DEBUG) << resultValues;

  // check that entries have expected values
  ASSERT_EQ(resultValues[0], 1.0);
  ASSERT_EQ(resultValues[1], 0.0);
  ASSERT_EQ(resultValues[48], 1.0);
  ASSERT_EQ(resultValues[98], 1.0);
  ASSERT_EQ(resultValues[146], 1.0);
}

TEST(FieldVariableTest, StructuredDeformable)
{
  // explicit functionSpace with node positions
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
  
  typedef FunctionSpace::FunctionSpace<MeshType,BasisFunctionType> FunctionSpaceType;
  typedef std::shared_ptr<FieldVariable::FieldVariableBaseFunctionSpace<FunctionSpaceType>> FieldVariableBaseFunctionSpaceType;
  typedef std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> FieldVariable1Type;
  typedef std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,2>> FieldVariable2Type;
  typedef std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> FieldVariable3Type;
  
  std::shared_ptr<FunctionSpaceType> functionSpace = finiteElementMethod.functionSpace();
  functionSpace->initialize();
  
  // 5x5 nodes, 2 dofs/node, 18 dofs/element, 50 dofs
  
  FieldVariableBaseFunctionSpaceType aBase = functionSpace->createFieldVariable("a", {"x","y"});
  FieldVariableBaseFunctionSpaceType bBase = functionSpace->createFieldVariable("b", 2);
  FieldVariableBaseFunctionSpaceType cBase = functionSpace->createFieldVariable("c");
  FieldVariable3Type d = functionSpace->template createFieldVariable<3>("d");
  
  
  FieldVariable2Type a = std::static_pointer_cast<FieldVariable::FieldVariable<FunctionSpaceType,2>>(aBase);
  FieldVariable2Type b = std::static_pointer_cast<FieldVariable::FieldVariable<FunctionSpaceType,2>>(bBase);
  FieldVariable1Type c = std::static_pointer_cast<FieldVariable::FieldVariable<FunctionSpaceType,1>>(cBase);
  
  // set all to 0.0
  a->setValues(0.0);
  
  ASSERT_EQ(a->getValue(0, 1), 0.0);
  ASSERT_EQ(a->getValue(0, 5), 0.0);
  ASSERT_EQ(a->getValue(1, 8), 0.0);
  
  // set all components, 1 dof
  Vec2 v0 = {1.0, 2.0};
  a->setValue(5, v0);
  
  ASSERT_EQ(a->getValue(0, 5), v0[0]);
  ASSERT_EQ(a->getValue(1, 5), v0[1]);
  
  // set all components, multiple dofs
  std::vector<Vec2> v1 = {Vec2({1.0, 2.0}), Vec2({3.0, 4.0}), Vec2({5.0, 6.0}), Vec2({7.0, 8.0}), Vec2({9.0, 10.0})};
  std::vector<dof_no_t> dofs = {4,8,2,7,10};

  a->setValues(dofs, v1);

  
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
  
  // set single component, 1 dof
  c->setValue(5, 5.0);
  ASSERT_EQ(c->getValue(0,5), 5.0);
  ASSERT_EQ(c->getValue(5), 5.0);
  ASSERT_EQ(c->getValue(0), 0.0);
  
  // set single component, multiple dofs
  std::vector<double> v2 = {1.0, 2.0, 3.0, 4.0, 5.0};
  c->setValues(dofs, v2);
  
  ASSERT_EQ(c->getValue(4), 1.0);
  ASSERT_EQ(c->getValue(8), 2.0);
  ASSERT_EQ(c->getValue(2), 3.0);
  ASSERT_EQ(c->getValue(7), 4.0);
  ASSERT_EQ(c->getValue(10), 5.0);
  
  std::vector<double> values10 = {11., 12., 13., 14.};
  std::vector<dof_no_t> dofs10 = {11, 12, 13, 14};
  c->setValues(dofs10, values10);
  
  /* values of c
     0.0, 0.0, 3.0, 0.0, 1.0,
     5.0, 0.0, 4.0, 2.0, 0.0,
     5.0,11.0,12.0,13.0,14.0,
     0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,
   */
// dofs element 1: 2,3,4, 7,8,9, 12,13,14
  
  std::array<double,FunctionSpaceType::nDofsPerElement()> values8;
  c->getElementValues(1.0, values8);
  
  std::array<double,FunctionSpaceType::nDofsPerElement()> reference8 = {
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
  a->getValuesWithoutGhosts(0, values);
  
  std::vector<double> reference0 = {
    0.0, 0.0, 5.0, 0.0, 1.0,
    1.0, 0.0, 7.0, 3.0, 0.0,
    9.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0
  };
  ASSERT_EQ(values, reference0);
  
  values.clear();
  a->getValuesWithoutGhosts(1, values);
  
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
  
  std::cout<<"--1";
  
  // all components, multiple dofs
  std::array<Vec2,3> values4;
  a->template getValues<3>(dofs2, values4);
  std::array<Vec2,3> reference4 = {Vec2({0.0,0.0}),Vec2({7.0,8.0}),Vec2({5.0,6.0})};
  ASSERT_EQ(values4,reference4);
  
  std::cout<<"--2";
  
//  node global nos.
//  10__11_12__13_14
//  ++==+==++==+==++
// 5||_6|_7||_8|__||9
//  ||__|__||__|__||
//   0  1  2   3  4
// dofs element 1: 2,3,4, 7,8,9, 12,13,14
  
  //! for a specific component, get the values corresponding to all element-local dofs
  // component 0
  std::array<double,FunctionSpaceType::nDofsPerElement()> values5;
  a->getElementValues(0, 1, values5);
  std::array<double,FunctionSpaceType::nDofsPerElement()> reference5 = {5.0, 0.0, 1.0, 7.0, 3.0, 0.0, 0.0, 0.0, 0.0};
  ASSERT_EQ(values5,reference5);
  
  std::cout<<"-3";
  // component 1
  std::array<double,FunctionSpaceType::nDofsPerElement()> values6;
  a->getElementValues(1, 1, values6);
  std::array<double,FunctionSpaceType::nDofsPerElement()> reference6 = {6.0, 0.0, 2.0, 8.0, 4.0, 0.0, 0.0, 0.0, 0.0};
  ASSERT_EQ(values6,reference6);
  
  
  
  std::cout<<"--4";
  //! get the values corresponding to all element-local dofs for all components
  std::array<Vec2,FunctionSpaceType::nDofsPerElement()> values7;
  a->getElementValues(1, values7);
  std::array<Vec2,FunctionSpaceType::nDofsPerElement()> reference7 = {Vec2({5.0, 6.0}), Vec2({0.0, 0.0}), Vec2({1.0, 2.0}), Vec2({7.0, 8.0}), Vec2({3.0, 4.0}), Vec2({0.0, 0.0}), Vec2({0.0, 0.0}), Vec2({0.0, 0.0}), Vec2({0.0, 0.0})};
  
  std::cout<<"--5";
  ASSERT_EQ(values7,reference7);
  
  // copy a to b
  b->setValues(*a);
  
  std::cout<<"--6";
  b->finishGhostManipulation();
  
  std::cout<<"--7";
  b->getElementValues(1, values7);
  ASSERT_EQ(values7,reference7);
  
  std::cout<<"--8";
  // set all values of a to 9.0
  a->setValues(9.0);
  
  std::cout<<"--9";
  std::array<Vec2,3> values9;
  a->getValues<3>(std::array<dof_no_t,3>{1,15,24}, values9);
  std::array<Vec2,3> reference9 = {Vec2({9.0, 9.0}), Vec2({9.0, 9.0}), Vec2({9.0, 9.0})};
  
  std::cout<<"--10";
  ASSERT_EQ(values9,reference9);
  /*
   * 
  //! copy the values from another field variable of the same type
  void setValues(FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents> &rhs);
  */
  
  /*
  //! create a non-geometry field field variable with no values being set, with given component names
  std::shared_ptr<FieldVariable::FieldVariableBaseFunctionSpace<FunctionSpace<MeshType,BasisFunctionType>>> createFieldVariable(std::string name, std::vector<std::string> componentNames);
  
  //! create a non-geometry field field variable with no values being set, with given number of components, the component names will be the numbers
  std::shared_ptr<FieldVariable::FieldVariableBaseFunctionSpace<FunctionSpace<MeshType,BasisFunctionType>>> createFieldVariable(std::string name, int nComponents=1);
  
  //! create a non-geometry field field variable with no values being set, with given number of components, the component names will be the numbers
  template <int nComponents>
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace<MeshType,BasisFunctionType>,nComponents>> createFieldVariable(std::string name);
  */
  
  /*
   * 
  //! set values for all components for dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
  void setValues(std::vector<dof_no_t> &dofGlobalNos, std::vector<std::array<double,nComponents>> &values, InsertMode petscInsertMode=INSERT_VALUES)
  
  //! set a single dof (all components) , after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
  void setValue(dof_no_t dofGlobalNo, std::array<double,nComponents> &value, InsertMode petscInsertMode=INSERT_VALUES)
 
   * 
  //! for a specific component, get all values
  //! @param onlyNodalValues: if this is true, for Hermite only the non-derivative values are retrieved
  void getValuesWithoutGhosts(int componentNo, std::vector<double> &values, bool onlyNodalValues=false);
  
  //! for a specific component, get values from their global dof no.s, as array, therefore templated by the number of elements, N, to retrieve
  template<int N>
  void getValuesWithoutGhosts(int componentNo, std::array<dof_no_t,N> dofGlobalNo, std::array<double,N> &values);
  
  //! for a specific component, get values from their global dof no.s, as vector
  void getValuesWithoutGhosts(int componentNo, std::vector<dof_no_t> dofGlobalNo, std::vector<double> &values);
  
  //! get values from their global dof no.s for all components
  template<int N>
  void getValuesWithoutGhosts(std::array<dof_no_t,N> dofGlobalNo, std::array<std::array<double,nComponents>,N> &values)
  
  //! for a specific component, get the values corresponding to all element-local dofs
  template<int N>
  void getElementValues(int componentNo, element_no_t elementNo, std::array<double,FunctionSpaceType::nDofsPerElement()> &values);
  
  //! get the values corresponding to all element-local dofs for all components
  void getElementValues(element_no_t elementNo, std::array<std::array<double,nComponents>,FunctionSpaceType::nDofsPerElement()> &values)
 
  //! for a specific component, get a single value from global dof no.
  double getValue(int componentNo, node_no_t dofGlobalNo);




  std::array<double,nComponents> getValue(node_no_t dofGlobalNo);
  
  
  
  
  
  //! get the values corresponding to all element-local dofs for all components
  void getElementValues(element_no_t elementNo, std::array<double,FunctionSpaceType::nDofsPerElement()> &values)
  
  //! get a single value from global dof no. for all components
  double getValue(node_no_t dofGlobalNo);
  
  //! set a single dof (all components) , after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
  void setValue(dof_no_t dofGlobalNo, double value, InsertMode petscInsertMode=INSERT_VALUES)

  //! set values for all components for dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
  void setValues(std::vector<dof_no_t> &dofGlobalNos, std::vector<double> &values, InsertMode petscInsertMode=INSERT_VALUES)
*/
}

TEST(FieldVariableTest, StructuredRegularFixed)
{
  // explicit functionSpace with node positions
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
  
  typedef FunctionSpace::FunctionSpace<MeshType,BasisFunctionType> FunctionSpaceType;
  typedef std::shared_ptr<FieldVariable::FieldVariableBaseFunctionSpace<FunctionSpaceType>> FieldVariableBaseFunctionSpaceType;
  typedef std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> FieldVariable1Type;
  typedef std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,2>> FieldVariable2Type;
  typedef std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> FieldVariable3Type;
  
  std::shared_ptr<FunctionSpaceType> functionSpace = std::static_pointer_cast<FunctionSpaceType>(finiteElementMethod.functionSpace());
  
  // 5x5 nodes, 2 dofs/node, 18 dofs/element, 50 dofs
  
  FieldVariableBaseFunctionSpaceType aBase = functionSpace->createFieldVariable("a", {"x","y"});
  FieldVariableBaseFunctionSpaceType bBase = functionSpace->createFieldVariable("b", 2);
  FieldVariableBaseFunctionSpaceType cBase = functionSpace->createFieldVariable("c");
  FieldVariable3Type d = functionSpace->template createFieldVariable<3>("d");
  
  
  FieldVariable2Type a = std::static_pointer_cast<FieldVariable::FieldVariable<FunctionSpaceType,2>>(aBase);
  FieldVariable2Type b = std::static_pointer_cast<FieldVariable::FieldVariable<FunctionSpaceType,2>>(bBase);
  FieldVariable1Type c = std::static_pointer_cast<FieldVariable::FieldVariable<FunctionSpaceType,1>>(cBase);
  
  // set all to 0.0
  a->setValues(0.0);
  
  ASSERT_EQ(a->getValue(0, 1), 0.0);
  ASSERT_EQ(a->getValue(0, 5), 0.0);
  ASSERT_EQ(a->getValue(1, 8), 0.0);
  
  // set all components, 1 dof
  Vec2 v0 = {1.0, 2.0};
  a->setValue(5, v0);
  
  ASSERT_EQ(a->getValue(0, 5), v0[0]);
  ASSERT_EQ(a->getValue(1, 5), v0[1]);
  
  // set all components, multiple dofs
  std::vector<Vec2> v1 = {Vec2({1.0, 2.0}), Vec2({3.0, 4.0}), Vec2({5.0, 6.0}), Vec2({7.0, 8.0}), Vec2({9.0, 10.0})};
  std::vector<dof_no_t> dofs = {4,8,2,7,10};
  a->setValues(dofs, v1);
  
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
  
  // set single component, 1 dof
  c->setValue(5, 5.0);
  ASSERT_EQ(c->getValue(0,5), 5.0);
  ASSERT_EQ(c->getValue(5), 5.0);
  ASSERT_EQ(c->getValue(0), 0.0);
  
  // set single component, multiple dofs
  std::vector<double> v2 = {1.0, 2.0, 3.0, 4.0, 5.0};
  c->setValues(dofs, v2);
  
  ASSERT_EQ(c->getValue(4), 1.0);
  ASSERT_EQ(c->getValue(8), 2.0);
  ASSERT_EQ(c->getValue(2), 3.0);
  ASSERT_EQ(c->getValue(7), 4.0);
  ASSERT_EQ(c->getValue(10), 5.0);
  
  std::vector<double> values10 = {11., 12., 13., 14.};
  std::vector<dof_no_t> dofs10 = {11, 12, 13, 14};
  c->setValues(dofs10, values10);
  
  /* values of c
     0.0, 0.0, 3.0, 0.0, 1.0,
     5.0, 0.0, 4.0, 2.0, 0.0,
     5.0,11.0,12.0,13.0,14.0,
     0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,
   */
// dofs element 1: 2,3,4, 7,8,9, 12,13,14
  
  std::array<double,FunctionSpaceType::nDofsPerElement()> values8;
  c->getElementValues(1.0, values8);
  
  std::array<double,FunctionSpaceType::nDofsPerElement()> reference8 = {
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
  a->getValuesWithoutGhosts(0, values);
  
  std::vector<double> reference0 = {
    0.0, 0.0, 5.0, 0.0, 1.0,
    1.0, 0.0, 7.0, 3.0, 0.0,
    9.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0
  };
  ASSERT_EQ(values, reference0);
  
  values.clear();
  a->getValuesWithoutGhosts(1, values);
  
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
  std::array<double,FunctionSpaceType::nDofsPerElement()> values5;
  a->getElementValues(0, 1, values5);
  std::array<double,FunctionSpaceType::nDofsPerElement()> reference5 = {5.0, 0.0, 1.0, 7.0, 3.0, 0.0, 0.0, 0.0, 0.0};
  ASSERT_EQ(values5,reference5);
  
  // component 1
  std::array<double,FunctionSpaceType::nDofsPerElement()> values6;
  a->getElementValues(1, 1, values6);
  std::array<double,FunctionSpaceType::nDofsPerElement()> reference6 = {6.0, 0.0, 2.0, 8.0, 4.0, 0.0, 0.0, 0.0, 0.0};
  ASSERT_EQ(values6,reference6);
  
  
  
  //! get the values corresponding to all element-local dofs for all components
  std::array<Vec2,FunctionSpaceType::nDofsPerElement()> values7;
  a->getElementValues(1, values7);
  std::array<Vec2,FunctionSpaceType::nDofsPerElement()> reference7 = {Vec2({5.0, 6.0}), Vec2({0.0, 0.0}), Vec2({1.0, 2.0}), Vec2({7.0, 8.0}), Vec2({3.0, 4.0}), Vec2({0.0, 0.0}), Vec2({0.0, 0.0}), Vec2({0.0, 0.0}), Vec2({0.0, 0.0})};
  
  ASSERT_EQ(values7,reference7);
  
  // copy a to b
  b->setValues(*a);
  b->getElementValues(1, values7);
  ASSERT_EQ(values7,reference7);
  
}

TEST(FieldVariableTest, UnstructuredDeformable)
{
  // explicit functionSpace with node positions
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
  
  typedef FunctionSpace::FunctionSpace<MeshType,BasisFunctionType> FunctionSpaceType;
  typedef std::shared_ptr<FieldVariable::FieldVariableBaseFunctionSpace<FunctionSpaceType>> FieldVariableBaseFunctionSpaceType;
  typedef std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> FieldVariable1Type;
  typedef std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,2>> FieldVariable2Type;
  typedef std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> FieldVariable3Type;
  
  std::shared_ptr<FunctionSpaceType> functionSpace = std::static_pointer_cast<FunctionSpaceType>(finiteElementMethod.functionSpace());
  
  // 5x5 nodes, 2 dofs/node, 18 dofs/element, 50 dofs
  
  FieldVariableBaseFunctionSpaceType aBase = functionSpace->createFieldVariable("a", {"x","y"});
  FieldVariableBaseFunctionSpaceType bBase = functionSpace->createFieldVariable("b", 2);
  FieldVariableBaseFunctionSpaceType cBase = functionSpace->createFieldVariable("c");
  FieldVariable3Type d = functionSpace->template createFieldVariable<3>("d");
  
  
  FieldVariable2Type a = std::static_pointer_cast<FieldVariable::FieldVariable<FunctionSpaceType,2>>(aBase);
  FieldVariable2Type b = std::static_pointer_cast<FieldVariable::FieldVariable<FunctionSpaceType,2>>(bBase);
  FieldVariable1Type c = std::static_pointer_cast<FieldVariable::FieldVariable<FunctionSpaceType,1>>(cBase);
  
  // set all to 0.0
  a->setValues(0.0);
  
  ASSERT_EQ(a->getValue(0, 1), 0.0);
  ASSERT_EQ(a->getValue(0, 5), 0.0);
  ASSERT_EQ(a->getValue(1, 8), 0.0);
  
  // set all components, 1 dof
  Vec2 v0 = {1.0, 2.0};
  a->setValue(5, v0);
  
  ASSERT_EQ(a->getValue(0, 5), v0[0]);
  ASSERT_EQ(a->getValue(1, 5), v0[1]);
  
  // set all components, multiple dofs
  std::vector<Vec2> v1 = {Vec2({1.0, 2.0}), Vec2({3.0, 4.0}), Vec2({5.0, 6.0}), Vec2({7.0, 8.0}), Vec2({9.0, 10.0})};
  std::vector<dof_no_t> dofs = {4,8,2,7,10};
  a->setValues(dofs, v1);
  
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
  
  // set single component, 1 dof
  c->setValue(5, 5.0);
  ASSERT_EQ(c->getValue(0,5), 5.0);
  ASSERT_EQ(c->getValue(5), 5.0);
  ASSERT_EQ(c->getValue(0), 0.0);
  
  // set single component, multiple dofs
  std::vector<double> v2 = {1.0, 2.0, 3.0, 4.0, 5.0};
  c->setValues(dofs, v2);
  
  ASSERT_EQ(c->getValue(4), 1.0);
  ASSERT_EQ(c->getValue(8), 2.0);
  ASSERT_EQ(c->getValue(2), 3.0);
  ASSERT_EQ(c->getValue(7), 4.0);
  ASSERT_EQ(c->getValue(10), 5.0);
  
  std::vector<double> values10 = {11., 12., 13., 14.};
  std::vector<dof_no_t> dofs10 = {11, 12, 13, 14};
  c->setValues(dofs10, values10);
  
  /* values of c
     0.0, 0.0, 3.0, 0.0, 1.0,
     5.0, 0.0, 4.0, 2.0, 0.0,
     5.0,11.0,12.0,13.0,14.0,
     0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,
   */
// dofs element 1: 2,9,10, 5,11,12, 8,13,14
  
  std::array<double,FunctionSpaceType::nDofsPerElement()> values8;
  c->getElementValues(1.0, values8);
  
  std::array<double,FunctionSpaceType::nDofsPerElement()> reference8 = {
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
  a->getValuesWithoutGhosts(0, values);
  
  std::vector<double> reference0 = {
    0.0, 0.0, 5.0, 0.0, 1.0,
    1.0, 0.0, 7.0, 3.0, 0.0,
    9.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0
  };
  ASSERT_EQ(values, reference0);
  
  values.clear();
  a->getValuesWithoutGhosts(1, values);
  
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
  std::array<double,FunctionSpaceType::nDofsPerElement()> values5;
  a->getElementValues(0, 1, values5);
  std::array<double,FunctionSpaceType::nDofsPerElement()> reference5 = {5.0, 0.0, 9.0, 1.0, 0.0, 0.0, 3.0, 0.0, 0.0};
  ASSERT_EQ(values5,reference5);
  
  // component 1
  std::array<double,FunctionSpaceType::nDofsPerElement()> values6;
  a->getElementValues(1, 1, values6);
  std::array<double,FunctionSpaceType::nDofsPerElement()> reference6 = {6.0, 0.0,10.0, 2.0, 0.0, 0.0, 4.0, 0.0, 0.0};
  ASSERT_EQ(values6,reference6);
  
  
  
  //! get the values corresponding to all element-local dofs for all components
  std::array<Vec2,FunctionSpaceType::nDofsPerElement()> values7;
  a->getElementValues(1, values7);
  std::array<Vec2,FunctionSpaceType::nDofsPerElement()> reference7 
    = {Vec2({5.0, 6.0}), Vec2({0.0, 0.0}), Vec2({9.0, 10.0}), Vec2({1.0, 2.0}), Vec2({0.0, 0.0}), Vec2({0.0, 0.0}), 
       Vec2({3.0, 4.0}), Vec2({0.0, 0.0}), Vec2({0.0, 0.0})};
  
  ASSERT_EQ(values7,reference7);
  
  // copy a to b
  b->setValues(*a);
  b->getElementValues(1, values7);
  ASSERT_EQ(values7,reference7);
  
}

}  // namespace
