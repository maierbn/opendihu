#include "data_management/solution_vector_mapping.h"

#include <vector>
#include "easylogging++.h"

SolutionVectorMapping::SolutionVectorMapping(bool canProvideInternalContiguousSolutionPointer) :
  scalingFactor_(1.0), canProvideInternalContiguousSolutionPointer_(canProvideInternalContiguousSolutionPointer)
{
}

void SolutionVectorMapping::setOutputRange(int outputIndexBegin, int outputIndexEnd)
{
  outputIndexBegin_ = outputIndexBegin;
  outputIndexEnd_ = outputIndexEnd;
  outputSize_ = outputIndexEnd_ - outputIndexBegin_;
}

void SolutionVectorMapping::setScalingFactor(double factor)
{
  scalingFactor_ = factor;
}


void SolutionVectorMapping::transfer(Vec &solution1, SolutionVectorMapping& solutionVectorMapping2, Vec &solution2)
{
  LOG(DEBUG) << "solution vector mapping (1): ["<<outputIndexBegin_<<","<<outputIndexEnd_<<"] (2): "
    <<"["<< solutionVectorMapping2.outputIndexBegin_<<","<<solutionVectorMapping2.outputIndexEnd_<<"]";
  if(canProvideInternalContiguousSolutionPointer_
    && solutionVectorMapping2.canProvideInternalContiguousSolutionPointer_)
  {
    LOG(DEBUG) << "transfer solution, both canProvideInternalContiguousSolutionPointer";
    
    const double *data1;
    VecGetArrayRead(solution1, &data1);
    
    LOG(DEBUG) << "data to be transferred ("<<outputSize_<<" entries): ";
    for(int i=outputIndexBegin_; i<outputIndexEnd_; i++)
      LOG(DEBUG) << "   " << data1[i];
    
    double *data2;
    VecGetArray(solution2, &data2);
    
    // copy data from data1 to data2
    memcpy(data2+solutionVectorMapping2.outputIndexBegin_, data1+outputIndexBegin_, outputSize_*sizeof(double));
    
    // scale data with scalingFactor
    if (scalingFactor_ != 1.0)
    {
      for(int i=0; i<outputSize_; i++)
      {
        data2[solutionVectorMapping2.outputIndexBegin_+i] *= scalingFactor_;
      }
    }
    
    VecRestoreArrayRead(solution1, &data1);
    VecRestoreArray(solution2, &data2);
  }
  else if (canProvideInternalContiguousSolutionPointer_ 
    && !solutionVectorMapping2.canProvideInternalContiguousSolutionPointer_)
  {
    LOG(DEBUG) << "transfer solution, only first canProvideInternalContiguousSolutionPointer";
    
    const double *data1;
    VecGetArrayRead(solution1, &data1);
    
    // copy data from data1 to solution2
    solutionVectorMapping2.copyFromDataPointerToSolution(data1, solution2);
    
    VecRestoreArrayRead(solution1, &data1);
  }
  else if(!canProvideInternalContiguousSolutionPointer_ 
    && solutionVectorMapping2.canProvideInternalContiguousSolutionPointer_)
  {
    LOG(DEBUG) << "transfer solution, only second canProvideInternalContiguousSolutionPointer";
    
    double *data2;
    VecGetArray(solution2, &data2);
    
    // copy data from solution1 to data2
    copyFromSolutionToDataPointer(solution1, data2);
    
    VecRestoreArray(solution2, &data2);
  }
  else if(!canProvideInternalContiguousSolutionPointer_ 
    && !solutionVectorMapping2.canProvideInternalContiguousSolutionPointer_)
  {
    LOG(DEBUG) << "transfer solution, none canProvideInternalContiguousSolutionPointer";
    
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
  
  // do scaling
}

void SolutionVectorMapping::copyFromSolutionToDataPointer(Vec& solution, double* data)
{
  // not yet implemented, because not yet needed.
  // This should handle case with strided access (for that a stride needs to be stored)
  
  // do scaling
}


