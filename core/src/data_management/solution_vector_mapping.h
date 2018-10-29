#pragma once

#include <Python.h>  // has to be the first included header
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
  SolutionVectorMapping();

  //! set the component no
  void setOutputComponentNo(int outputComponentNo);

  //! set the index range of the solution vector to be considered for transfor
  //void setOutputRange(int outputIndexBegin, int outputIndexEnd);

  //! set a scaling factor by which the result will be scaled when transferring
  void setScalingFactor(double factor);

  //! transfer data from solution1 to solution2, where solution1 corresponds to "this" object, solution2 corresponds to solutionVectorMapping2
  template<typename FieldVariable1, typename FieldVariable2>
  void transfer(FieldVariable1 &solution1, std::shared_ptr<SolutionVectorMapping> solutionVectorMapping2, FieldVariable2 &solution2);

  //void copyFromSolutionToDataPointer(Vec &solution, double *data);
  //void copyFromDataPointerToSolution(const double *data, Vec &solution);
private:

 //int outputIndexBegin_;
 //int outputIndexEnd_;
 //int outputSize_;
 double scalingFactor_;   ///< the prefactor with which the results is scaled when transferring from solution1 to solution2

 int outputComponentNo_;  ///< the component of the solution vector which should be transferred to the other solution vector

 //bool canProvideInternalContiguousSolutionPointer_;
};

#include "solution_vector_mapping.tpp"
