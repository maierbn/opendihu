#pragma once

#include <petscmat.h>
#include <memory>

#include "data_management/data.h"
#include "control/types.h"
#include "mesh/mesh.h"

/**
 * The Data classes contain each a vector that stores the solution. Often, the values need to be accessed to 
 * continue computation using them, e.g. in operator splittings. Then possibly not all the values need to be accessed
 * but only some of them while the others are only needed for the own computation. An example is a cellml model that 
 * contains lots of internal states and only a single state variable is needed in the diffusion equation.
 * This class stores information for a single solution vector/equation object on how the needed values can be extracted.
 */
class SolutionVectorMapping
{
public:
  //! constructor
  SolutionVectorMapping(bool canProvideInternalContiguousSolutionPointer);
 
  //! set the index range of the solution vector to be considered for transfor
  void setOutputRange(int outputIndexBegin, int outputIndexEnd);
  
  //! set a scaling factor by which the result will be scaled when transferring
  void setScalingFactor(double factor);
  
  //! transfer data from solution1 to solution2, where solution1 corresponds to "this" object, solution2 corresponds to solutionVectorMapping2
  void transfer(Vec &solution1, SolutionVectorMapping &solutionVectorMapping2, Vec &solution2);
  
  void copyFromSolutionToDataPointer(Vec &solution, double *data);
  void copyFromDataPointerToSolution(const double *data, Vec &solution);
private:
 
 int outputIndexBegin_;
 int outputIndexEnd_;
 int outputSize_;
 double scalingFactor_;   ///< the prefactor with which the results is scaled when transferring from solution1 to solution2
 
 bool canProvideInternalContiguousSolutionPointer_;
};