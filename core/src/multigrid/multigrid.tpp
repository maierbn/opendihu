#include "multigrid/multigrid.h"
#include "../field_variable/structured/04_field_variable_set_get_component_dependent_structured.h"

#include "utility/python_utility.h"
#include "control/performance_measurement.h"

#include "function_space/function_space.h"
#include "field_variable/field_variable.h"
#include "control/types.h"

namespace Multigrid
{

template<typename FiniteElement1, typename FiniteElement2>
Multigrid<FiniteElement1, FiniteElement2>::Multigrid(DihuContext context) :
finiteElement1_(context["Term1"]), finiteElement2_(context["Term2"]), initialized_(false)
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
  solution_1=finiteElement1_.data().solution();  
  //std::array<double,FunctionSpace1::nDofsPerElement()> values_1;
  //element_no_t elementNo_1;
  //const std::vector<dof_no_t> dofLocalNo_1;
  std::vector<double> values_1;
  
  typedef typename FiniteElement2::FunctionSpace FunctionSpace2;
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace2,1>>
  solution_2=finiteElement2_.data().solution();
  //std::array<double,FunctionSpace2::nDofsPerElement()> values_2;
  //element_no_t elementNo_2;
  std::vector<double> values_2;
  
  
	for (int i = 0; i<numCycles_; i++)
	{
		finiteElement1_.solveMG();
		
    //transfer
    //--------
    // in case we need access to each element
    //FieldVariable::FieldVariableSetGetComponent<FunctionSpace1,1>::
    //getElementValues(elementNo_1, &values_1); 
    
    solution_1->getValuesWithoutGhosts(values_1);
    
    solution_2->setValuesWithoutGhosts(values_1);
    solution_2->startGhostManipulation();    
    solution_2->finishGhostManipulation();
    
    finiteElement2_.solveMG();
		//transfer
    //--------
    solution_2->getValuesWithoutGhosts(values_2);
    
    solution_1->setValuesWithoutGhosts(values_2);
    solution_1->startGhostManipulation();    
    solution_1->finishGhostManipulation();
    
		}
}		

};    // namespace
