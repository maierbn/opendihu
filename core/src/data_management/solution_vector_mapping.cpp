#include "data_management/solution_vector_mapping.h"

#include <vector>

SolutionVectorMapping::SolutionVectorMapping(bool canProvideInternalContiguousSolutionPointer) :
  canProvideInternalContiguousSolutionPointer_(canProvideInternalContiguousSolutionPointer)
{
}

void SolutionVectorMapping::setOutputRange(int outputIndexBegin, int outputIndexEnd)
{
  outputIndexBegin_ = outputIndexBegin;
  outputIndexEnd_ = outputIndexEnd;
  outputSize_ = outputIndexEnd_ - outputIndexBegin_;
}

void SolutionVectorMapping::transfer(Vec &solution1, SolutionVectorMapping& solutionVectorMapping2, Vec &solution2)
{
  if(canProvideInternalContiguousSolutionPointer_
    && solutionVectorMapping2.canProvideInternalContiguousSolutionPointer_)
  {
    const double *data1;
    VecGetArrayRead(solution1, &data1);
    
    double *data2;
    VecGetArray(solution2, &data2);
    
    // copy data from data1 to data2
    memcpy(data2+solutionVectorMapping2.outputIndexBegin_, data1+outputIndexBegin_, outputSize_*sizeof(double));
    
    VecRestoreArrayRead(solution1, &data1);
    VecRestoreArray(solution2, &data2);
  }
  else if (canProvideInternalContiguousSolutionPointer_ 
    && !solutionVectorMapping2.canProvideInternalContiguousSolutionPointer_)
  {
    const double *data1;
    VecGetArrayRead(solution1, &data1);
    
    // copy data from data1 to solution2
    solutionVectorMapping2.copyFromDataPointerToSolution(data1, solution2);
    
    VecRestoreArrayRead(solution1, &data1);
  }
  else if(!canProvideInternalContiguousSolutionPointer_ 
    && solutionVectorMapping2.canProvideInternalContiguousSolutionPointer_)
  {
    double *data2;
    VecGetArray(solution2, &data2);
    
    // copy data from solution1 to data2
    copyFromSolutionToDataPointer(solution1, data2);
    
    VecRestoreArray(solution2, &data2);
  }
  else if(!canProvideInternalContiguousSolutionPointer_ 
    && !solutionVectorMapping2.canProvideInternalContiguousSolutionPointer_)
  {
    std::vector<double> data;
    data.reserve(outputSize_);
    
    copyFromSolutionToDataPointer(solution1, data.data());
    solutionVectorMapping2.copyFromDataPointerToSolution(data.data(), solution2);
  }
}

void SolutionVectorMapping::copyFromDataPointerToSolution(const double* data, Vec& solution)
{
  // not yet implemented, because not yet needed.
  // This should handle case with strided access (for that a stride needs to be stored)
}

void SolutionVectorMapping::copyFromSolutionToDataPointer(Vec& solution, double* data)
{
  // not yet implemented, because not yet needed.
  // This should handle case with strided access (for that a stride needs to be stored)
}


