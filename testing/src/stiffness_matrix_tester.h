#pragma once
  
#include <vector>
#include <map>
#include "control/petsc_utility.h"

namespace SpatialDiscretization
{
  
class StiffnessMatrixTester
{
public:
  template<class MeshType, class BasisFunctionType, class EquationType>
  static void compareMatrix(
    FiniteElementMethod<MeshType, BasisFunctionType, EquationType> &finiteElementMethod,
    std::vector<double> &referenceMatrix
                           )
  {
    Mat &stiffnessMatrix = finiteElementMethod.data_.stiffnessMatrix();
    std::vector<double> matrix;
    PetscUtility::getMatrixEntries(stiffnessMatrix, matrix);
    
    ASSERT_EQ(matrix.size(), referenceMatrix.size()) << "Matrix has wrong number of entries";
    for(unsigned int i=0; i<matrix.size(); i++)
    {
      EXPECT_EQ(matrix[i], referenceMatrix[i]) << "Matrix entry no. " << i << " differs";
    }
  }
  
  template<class MeshType, class BasisFunctionType, class EquationType>
  static void compareRhs(
    FiniteElementMethod<MeshType, BasisFunctionType, EquationType> &finiteElementMethod,
    std::vector<double> &referenceRhs
                           )
  {
    Vec &rhsVector = finiteElementMethod.data_.rightHandSide();
    std::vector<double> rhs;
    PetscUtility::getVectorEntries(rhsVector, rhs);
    
    ASSERT_EQ(rhs.size(), referenceRhs.size()) << "Rhs has wrong number of entries";
    for(unsigned int i=0; i<rhs.size(); i++)
    {
      EXPECT_EQ(rhs[i], referenceRhs[i]) << "Rhs entry no. " << i << " differs";
    }
  }
  
  template<class MeshType, class BasisFunctionType, class EquationType>
  static void checkDirichletBCInSolution(
    FiniteElementMethod<MeshType, BasisFunctionType, EquationType> &finiteElementMethod,
    std::map<int, double> &dirichletBC)
  {
    Vec &solutionVector = finiteElementMethod.data_.solution();
    std::vector<double> solution;
    PetscUtility::getVectorEntries(solutionVector, solution);
    
    for(std::pair<int, double> entry : dirichletBC)
    {
      double difference = fabs(solution[entry.first]-entry.second);
      EXPECT_LE(difference, 1e-15) 
        << "Dirichlet BC on node " << entry.first << ", value "<<entry.second << " is not met! "
        << "Actual value: "<<solution[entry.first] << ", Difference: " << difference;
    }
  }
  
  template<class MeshType, class BasisFunctionType, class EquationType>
  static void checkEqual(
    FiniteElementMethod<MeshType, BasisFunctionType, EquationType> &finiteElementMethod1,
    FiniteElementMethod<MeshType, BasisFunctionType, EquationType> &finiteElementMethod2)
  {
    // rhsVector
    Vec &rhsVector1 = finiteElementMethod1.data_.rightHandSide();
    std::vector<double> rhs1;
    PetscUtility::getVectorEntries(rhsVector1, rhs1);
    Vec &rhsVector2 = finiteElementMethod2.data_.rightHandSide();
    std::vector<double> rhs2;
    PetscUtility::getVectorEntries(rhsVector2, rhs2);
    
    ASSERT_EQ(rhs1.size(), rhs2.size()) << "Rhs has wrong number of entries";
    for(unsigned int i=0; i<rhs1.size(); i++)
    {
      EXPECT_EQ(rhs1[i], rhs2[i]) << "Rhs entry no. " << i << " differs";
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
      EXPECT_EQ(matrix1[i], matrix2[i]) << "Matrix entry no. " << i << " differs";
    }
    
    
  }

  template<class MeshType, class BasisFunctionType>
  static void testRhsDiscretizationMatrix(
    FiniteElementMethod<MeshType, BasisFunctionType, Equation::Static::Poisson> &finiteElementMethod1,
    FiniteElementMethod<MeshType, BasisFunctionType, Equation::Dynamic::Diffusion> &finiteElementMethod2,
    std::vector<double> &rhsValues
  )
  {
    // create the discretization matrix if it does not already exist
    finiteElementMethod2.createRhsDiscretizationMatrix();
    Mat &dmatrix = finiteElementMethod2.data_.discretizationMatrix();
    
    int n, m;
    MatGetSize(dmatrix, &n, &m);
    std::cout << "matrix size: " << n << "x" << m << std::endl;
    Vec rhsStrong, rhsWeak;
      
    PetscUtility::createVector(rhsStrong, n, "rhs strong");
    PetscUtility::createVector(rhsWeak, n, "rhs weak");
    
    PetscUtility::setVector(rhsValues, rhsStrong);

    LOG(DEBUG) << "discretizationMatrix: " << PetscUtility::getStringMatrixVector(dmatrix, rhsStrong);    
    MatMult(dmatrix, rhsStrong, rhsWeak);
    
    LOG(DEBUG) << "dmatrix * rhsStrong = rhsWeak: " << PetscUtility::getStringVector(rhsWeak);
    
    // dmatrix * f_strong = rhs_weak
    Vec &rhs = finiteElementMethod1.data_.rightHandSide();   // rhs in weak formulation
   
    LOG(DEBUG) << "using stencil: "<<PetscUtility::getStringVector(rhs);
    
    std::vector<double> rhsWeakDMatrix, rhsWeakStencil;
    PetscUtility::getVectorEntries(rhsWeak, rhsWeakDMatrix);
    PetscUtility::getVectorEntries(rhs, rhsWeakStencil);
    
    // compare vectors
    ASSERT_EQ(rhsWeakDMatrix.size(), rhsWeakStencil.size());
    
    for(unsigned int i=0; i<rhsWeakDMatrix.size(); i++)
    {
      double difference = fabs(rhsWeakDMatrix[i]-rhsWeakStencil[i]);
      EXPECT_LE(difference, 1e-15) 
        << "Rhs entry number " << i << " is different. Using dmatrix: " << rhsWeakDMatrix[i] 
        << ", using stencil: " << rhsWeakStencil[i] << ", Difference: " << difference;
    }
  }
  
};

};