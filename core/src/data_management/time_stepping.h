#pragma once

#include <petscmat.h>
#include <memory>

#include "control/types.h"
#include "mesh/mesh.h"
#include "data_management/data.h"

class DihuContext;

namespace Data
{
 
class TimeStepping : public Data
{
public:
 
  //! constructor
  TimeStepping(const DihuContext &context);
  
  //! destructur
  ~TimeStepping();
  
  //! return a reference to the solution vector
  Vec &solution();
  
  //! return a reference to the increment vector
  Vec &increment();
 
  //! print all stored data to stdout
  void print();
  
private:
 
  //! initializes the vectors with size
  void createPetscObjects();
  
  bool disablePrinting_ = false;    ///< if printing vectors is disabled
  
  Vec solution_;            ///< the vector of the variable of interest
  Vec increment_;        ///< the vector for delta u
};

} // namespace Data