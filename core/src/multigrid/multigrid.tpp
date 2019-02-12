#include "multigrid/multigrid.h"
#include "../field_variable/structured/04_field_variable_set_get_component_dependent_structured.h"

#include "utility/python_utility.h"
#include "control/performance_measurement.h"

#include "function_space/function_space.h"
#include "field_variable/field_variable.h"
#include "control/types.h"
#include "multigrid/transfer.h"

namespace Multigrid
{

template<typename FiniteElement1, typename FiniteElement2>
Multigrid<FiniteElement1, FiniteElement2>::Multigrid(DihuContext context) :
finiteElement1_(context["Multigrid"]["Term1"]), finiteElement2_(context["Multigrid"]["Term2"]), initialized_(false),
numCycles_(100)
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
  typedef typename FiniteElement1::FunctionSpace FunctionSpace1;  
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace1,1>>
  solution1=finiteElement1_.data().solution();  
  //std::array<double,FunctionSpace1::nDofsPerElement()> values1;
  //element_no_t elementNo1;
  //const std::vector<dof_no_t> dofLocalNo1;
  std::vector<double> values1;
  
  typedef typename FiniteElement2::FunctionSpace FunctionSpace2;
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace2,1>>
  solution2=finiteElement2_.data().solution();
  //std::array<double,FunctionSpace2::nDofsPerElement()> values2;
  //element_no_t elementNo2;
  std::vector<double> values2;
  
  
 for (int i = 0; i<numCycles_; i++)
 {
  finiteElement1_.solveMG();

    //transfer
    //--------
    // in case we need access to each element
    //FieldVariable::FieldVariableSetGetComponent<FunctionSpace1,1>::
    //getElementValues(elementNo1, &values1);
    
    values1.clear();
    solution1->getValuesWithoutGhosts(values1);
    
    //restriction
    Transfer::restriction(&values1);
   
    solution2->setValuesWithoutGhosts(values1);
    //solution2->startGhostManipulation();    
    //solution2->finishGhostManipulation();
    
    finiteElement2_.solveMG();
    //transfer
    //--------
    values2.clear();
    solution2->getValuesWithoutGhosts(values2);
    
    //prolongation
    Transfer::prolongation(&values2);
    
    solution1->setValuesWithoutGhosts(values2);
    //solution1->startGhostManipulation();    
    //solution1->finishGhostManipulation();
    
  }
}

};    // namespace
