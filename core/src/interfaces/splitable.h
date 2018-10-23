#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>

#include "data_management/solution_vector_mapping.h"

/**
 *  Class that can be part of a splitting scheme
 */
class Splitable
{
public:
  //! constructor
  Splitable();

  //! return the solution vector mapping object, that contains information on if there are more internal values stored in the data_ object than may be needed for further computationo
  std::shared_ptr<SolutionVectorMapping> solutionVectorMapping();

protected:
  std::shared_ptr<SolutionVectorMapping> solutionVectorMapping_;   ///< the solution vector mapping object that contains information if for further computation only a subset of the stored entries in the data_.solution vector will be needed
};
