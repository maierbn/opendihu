#include <iostream>
#include <cstdlib>

#include <iostream>
#include <petscksp.h>
#include <vector>
#include <iomanip>
#include <numeric>

#include "control/dihu_context.h"
#include "control/simulation.h"
#include "control/types.h"
#include "easylogging++.h"

void print(Mat &stiffnessMatrix_)
{
 
  LOG(INFO)<<"======================";
  int nRows, nColumns;
  MatGetSize(stiffnessMatrix_, &nRows, &nColumns);
  
  LOG(INFO)<<"stiffnessMatrix ("<<nRows<<" x "<<nColumns<<"):";
  
  std::vector<int> rowIndices(nRows);
  std::iota(rowIndices.begin(), rowIndices.end(), 0);
  std::vector<int> columnIndices(nColumns);
  std::iota(columnIndices.begin(), columnIndices.end(), 0);
  std::vector<double> matrixValues(nRows*nColumns);
  
  MatGetValues(stiffnessMatrix_, nRows, rowIndices.data(), nColumns, columnIndices.data(), matrixValues.data());
  
  std::stringstream s;
  for (int i=0; i<nRows; i++)
  {
    s<<std::setw(3)<<std::setfill(' ')<<i<<" ";
    for (int j=0; j<nColumns; j++)
    {
      s<<matrixValues[i*nRows + j]<<" ";
    }
    s<<std::endl;
  }
  s<<std::endl;
  LOG(INFO) << std::endl<<s.str();
  
  MatInfo info;
  MatGetInfo(stiffnessMatrix_, MAT_LOCAL, &info);
  
  LOG(INFO)<<"Matrix info: "<<std::endl
    <<"block_size: "<<info.block_size<<std::endl
    <<"number of nonzeros: allocated: "<<info.nz_allocated<<", used: "<<info.nz_used<<", unneeded: "<<info.nz_unneeded
    <<"memory allocated: "<<info.memory<<std::endl
    <<"number of matrix assemblies called: "<<info.assemblies<<std::endl
    <<"number of mallocs during MatSetValues(): "<<info.mallocs<<std::endl
    <<"fill ratio for LU/ILU: given: "<<info.fill_ratio_given<<", needed: "<<info.fill_ratio_needed<<std::endl 
    <<"number of mallocs during factorization: "<<info.factor_mallocs<<std::endl;
    
    
  LOG(INFO)<<"======================";
   
}

void test()
{
  // testing how to use PETSc matrices
  Mat stiffnessMatrix_;
  Mat &stiffnessMatrix = stiffnessMatrix_;

  
  // PETSc MatCreateAIJ parameters
  const int d_nz = 3;   // number of nonzeros per row in DIAGONAL portion of local submatrix (same value is used for all local rows)
  const int o_nz = 0;   //  number of nonzeros per row in the OFF-DIAGONAL portion of local submatrix (same value is used for all local rows)
  
  int n = 5;
  int nElements = n;
  
  PetscErrorCode ierr;
  ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, n, n, 
                      d_nz, NULL, o_nz, NULL, &stiffnessMatrix_); CHKERRV(ierr);
  ierr = MatMPIAIJSetPreallocation(stiffnessMatrix_, d_nz, NULL, o_nz, NULL); CHKERRV(ierr);
  
  LOG(DEBUG) << "set values";
  double elementLength = 0.1;
  
  for (node_idx_t elementNo = 0; elementNo < nElements; elementNo++)
  {
    // stencil for -Δu in 1D: [1 -2 1]
   
    //                 matrix           row        column
    ierr = MatSetValue(stiffnessMatrix, elementNo, elementNo, -2.0*elementLength, INSERT_VALUES); CHKERRV(ierr);
   
    if (elementNo+1 < nElements)
    {
     ierr = MatSetValue(stiffnessMatrix, elementNo, elementNo+1, 1.0*elementLength, INSERT_VALUES); CHKERRV(ierr);
    }
    if (elementNo-1 >= 0)
    {
     ierr = MatSetValue(stiffnessMatrix, elementNo, elementNo-1, 1.0*elementLength, INSERT_VALUES); CHKERRV(ierr);
    }
  }
  
  LOG(DEBUG) << "MatAssemblyBegin";
  ierr = MatAssemblyBegin(stiffnessMatrix_, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(stiffnessMatrix_, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);

  print(stiffnessMatrix);

  node_idx_t boundaryConditionNodeIndex = 3;
  double boundaryConditionValue = 5.0;
  
  // get the column boundaryConditionNodeIndex of the stiffness matrix. It is needed for updating the rhs
  std::vector<int> rowIndices(nElements);
  std::iota (rowIndices.begin(), rowIndices.end(), 0);    // fill with increasing numbers: 0,1,2,...
  std::vector<int> columnIndices = {boundaryConditionNodeIndex};
  
  std::vector<double> coefficients(nElements);
  
  LOG(DEBUG) << "MatGetValues";
  ierr = MatGetValues(stiffnessMatrix, nElements, rowIndices.data(), 1, columnIndices.data(), coefficients.data());
  
  LOG(DEBUG) << "coefficients: ";
  for(int i=0; i<coefficients.size(); i++)
    LOG(DEBUG) << coefficients[i];
  
  LOG(DEBUG) << "MatZeroRowsColumns";
  MatZeroRowsColumns(stiffnessMatrix, 1, &boundaryConditionNodeIndex, 1.0, NULL, NULL);
  
  print(stiffnessMatrix);
}

int main(int argc, char *argv[])
{
  std::cout << "hu" << std::endl;
  
  DihuContext dihuContext(argc, argv);
  Simulation simulation;
  simulation.debug();
  
  test();
  return EXIT_SUCCESS;
}