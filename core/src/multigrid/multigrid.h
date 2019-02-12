#pragma once

#include <Python.h>  // has to be the first included header
#include "control/dihu_context.h"
#include "output_writer/manager.h"
#include "interfaces/runnable.h"

#include "easylogging++.h"


namespace Multigrid
{

template<typename FiniteElement1, typename FiniteElement2>
class Multigrid : public Runnable
{
public:

  //! constructor
  Multigrid(DihuContext context);
  
  //! run the simulation
  void run();
  
   //! return a solution vector
  Vec &solution();
  
  //! initialize data
  void initialize();
  
  //! reset state such that new initialization becomes necessary
  virtual void reset();

  //! return the data object
  //Data &data();
  
  void solveMG();
  
protected:

FiniteElement1 finiteElement1_;
FiniteElement2 finiteElement2_;

DihuContext context_;
PythonConfig specificSettings_;

bool initialized_;
int numCycles_ ;
std::string multigridType_;
};

}  // namespace

#include "multigrid/multigrid.tpp"
