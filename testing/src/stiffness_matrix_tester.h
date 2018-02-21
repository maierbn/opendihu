#pragma once
  
#include <Python.h>  // this has to be the first included header
#include <vector>
#include <map>
#include "utility/petsc_utility.h"

namespace SpatialDiscretization
{
  
class StiffnessMatrixTester
{
public:
  template<typename MeshType, typename BasisFunctionType, typename IntegratorType, typename EquationType>
  static void compareMatrix(
    FiniteElementMethod<MeshType, BasisFunctionType, IntegratorType, EquationType> &finiteElementMethod,
    std::vector<double> &referenceMatrix
                           )
  {
    Mat &stiffnessMatrix = finiteElementMethod.data_.stiffnessMatrix();
    std::vector<double> matrix;
    PetscUtility::getMatrixEntries(stiffnessMatrix, matrix);
    
    ASSERT_EQ(matrix.size(), referenceMatrix.size()) << "Matrix has wrong number of entries";
    for(unsigned int i=0; i<matrix.size(); i++)
    {
      double difference = fabs(matrix[i]-referenceMatrix[i]);
      EXPECT_LE(difference, 1e-14) << "Matrix entry no. " << i << " differs by " << difference << ".";
    }
  }
  
  template<typename MeshType, typename BasisFunctionType, typename IntegratorType, typename EquationType>
  static void compareRhs(
    FiniteElementMethod<MeshType, BasisFunctionType, IntegratorType, EquationType> &finiteElementMethod,
    std::vector<double> &referenceRhs
                           )
  {
    Vec &rhsVector = finiteElementMethod.data_.rightHandSide().values();
    std::vector<double> rhs;
    PetscUtility::getVectorEntries(rhsVector, rhs);
    
    ASSERT_EQ(rhs.size(), referenceRhs.size()) << "Rhs has wrong number of entries";
    for(unsigned int i=0; i<rhs.size(); i++)
    {
      double difference = fabs(rhs[i] - referenceRhs[i]);
      EXPECT_LE(difference, 1e-14) << "Rhs entry no. " << i << " differs by " << difference << ".";
    }
  }
  
  template<typename MeshType, typename BasisFunctionType, typename IntegratorType, typename EquationType>
  static void checkDirichletBCInSolution(
    FiniteElementMethod<MeshType, BasisFunctionType, IntegratorType, EquationType> &finiteElementMethod,
    std::map<int, double> &dirichletBC)
  {
    Vec &solutionVector = finiteElementMethod.data_.solution().values();
    std::vector<double> solution;
    PetscUtility::getVectorEntries(solutionVector, solution);
    
    for(std::pair<int, double> entry : dirichletBC)
    {
      double difference = fabs(solution[entry.first]-entry.second);
      EXPECT_LE(difference, 1e-14) 
        << "Dirichlet BC on node " << entry.first << ", value "<<entry.second << " is not met! "
        << "Actual value: "<<solution[entry.first] << ", Difference: " << difference;
    }
  }
  
  template<typename MeshType1, typename BasisFunctionType1, typename IntegratorType1, typename EquationType1,
           typename MeshType2, typename BasisFunctionType2, typename IntegratorType2, typename EquationType2>
  static void checkEqual(
    FiniteElementMethod<MeshType1, BasisFunctionType1, IntegratorType1, EquationType1> &finiteElementMethod1,
    FiniteElementMethod<MeshType2, BasisFunctionType2, IntegratorType2, EquationType2> &finiteElementMethod2)
  {
    // rhsVector
    Vec &rhsVector1 = finiteElementMethod1.data_.rightHandSide().values();
    std::vector<double> rhs1;
    PetscUtility::getVectorEntries(rhsVector1, rhs1);
    Vec &rhsVector2 = finiteElementMethod2.data_.rightHandSide().values();
    std::vector<double> rhs2;
    PetscUtility::getVectorEntries(rhsVector2, rhs2);
    
    ASSERT_EQ(rhs1.size(), rhs2.size()) << "Rhs has wrong number of entries";
    for(unsigned int i=0; i<rhs1.size(); i++)
    {
      double difference = fabs(rhs1[i]-rhs2[i]);
      EXPECT_LE(difference, 1e-14) << "Rhs entry no. " << i << " differs by " << difference << ".";
    }
    
    // stiffness matrix 
    Mat &stiffnessMatrix1 = finiteElementMethod1.data_.stiffnessMatrix();
    std::vector<double> matrix1;
    PetscUtility::getMatrixEntries(stiffnessMatrix1, matrix1);
    Mat &stiffnessMatrix2 = finiteElementMethod2.data_.stiffnessMatrix();
    std::vector<double> matrix2;
    PetscUtility::getMatrixEntries(stiffnessMatrix2, matrix2);
    
    ASSERT_EQ(matrix1.size(), matrix2.size()) << "Matrix has wrong number of entries";
    for(unsigned int i=0; i<matrix1.size(); i++)
    {
      double difference = fabs(matrix1[i]-matrix2[i]);
      EXPECT_LE(difference, 2e-14) << "Matrix entry no. " << i << " differs by " << difference << ".";
    }
    
    
  }

  template<typename T1, typename T2>
  static void testRhsDiscretizationMatrix(
    T1 &finiteElementMethod1,
    T2 &finiteElementMethod2,
    std::vector<double> &rhsValues
  )
  {
    // create the discretization matrix if it does not already exist
    finiteElementMethod2.setRhsDiscretizationMatrix();
    Mat &massMatrix = finiteElementMethod2.data_.massMatrix();
    
    int n, m;
    MatGetSize(massMatrix, &n, &m);
    LOG(DEBUG) << "matrix size: " << n << "x" << m << std::endl;
    Vec rhsStrong, rhsWeak;
      
    PetscUtility::createVector(rhsStrong, n, "rhs strong");
    PetscUtility::createVector(rhsWeak, n, "rhs weak");
    
    PetscUtility::setVector(rhsValues, rhsStrong);

    LOG(DEBUG) << "massMatrix: " << PetscUtility::getStringMatrixVector(massMatrix, rhsStrong);    
    MatMult(massMatrix, rhsStrong, rhsWeak);
    
    LOG(DEBUG) << "massMatrix * rhsStrong = rhsWeak: " << PetscUtility::getStringVector(rhsWeak);
    
    // massMatrix * f_strong = rhs_weak
    Vec &rhs = finiteElementMethod1.data_.rightHandSide().values();   // rhs in weak formulation
   
    LOG(DEBUG) << "using stencil: "<<PetscUtility::getStringVector(rhs);
    
    std::vector<double> rhsWeakDMatrix, rhsWeakStencil;
    PetscUtility::getVectorEntries(rhsWeak, rhsWeakDMatrix);
    PetscUtility::getVectorEntries(rhs, rhsWeakStencil);
    
    // compare vectors
    ASSERT_EQ(rhsWeakDMatrix.size(), rhsWeakStencil.size());
    
    for(unsigned int i=0; i<rhsWeakDMatrix.size(); i++)
    {
      double difference = fabs(rhsWeakDMatrix[i]-rhsWeakStencil[i]);
      EXPECT_LE(difference, 1e-14) 
        << "Rhs entry number " << i << " is different. Using massMatrix: " << rhsWeakDMatrix[i] 
        << ", using stencil: " << rhsWeakStencil[i] << ", Difference: " << difference;
    }
  }
  
};

};
