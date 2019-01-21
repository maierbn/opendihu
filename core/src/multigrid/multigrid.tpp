#include "multigrid/multigrid.h"

#include "utility/python_utility.h"
#include "control/performance_measurement.h"


namespace Multigrid
{

template<typename FiniteElement1, typename FiniteElement2>
Multigrid<FiniteElement1, FiniteElement2>::
Multigrid(DihuContext context)
  {
	
}
template<typename FiniteElement1, typename FiniteElement2>
void Multigrid<FiniteElement1, FiniteElement2>::
reset()
{
	initialized_ = false;
	finiteElement1_.reset();
	finiteElement2_.reset();
}
template<typename FiniteElement1, typename FiniteElement2>
void Multigrid<FiniteElement1, FiniteElement2>::
initialize()
{
  if (initialized_)
    return;
  LOG(TRACE) << "  Multigrid::initialize";
    // initialize Finite Element objects
    LOG(DEBUG) << "  Multigrid::initialize FiniteElement1";
    finiteElement1_.initialize();
    LOG(DEBUG) << "  Multigrid::initialize FiniteElement2";
    finiteElement2_.initialize();
  
    initialized_ = true;
}

template<typename FiniteElement1, typename FiniteElement2>
void Multigrid<FiniteElement1, FiniteElement2>::
run()
{
  // initialize data structurures
  this->initialize();
  //solveMG here or use FEM Solve for each Term?
  solveMG();

}
template<typename FiniteElement1, typename FiniteElement2>
void Multigrid<FiniteElement1, FiniteElement2>::
solveMG()
{
	
}

};    // namespace
