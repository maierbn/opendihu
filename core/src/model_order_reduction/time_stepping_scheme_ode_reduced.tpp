#include "model_order_reduction/time_stepping_scheme_ode_reduced.h"

#include <Python.h>
#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include<petscmat.h>
#include "mesh/mesh_manager/mesh_manager.h"
#include "function_space/function_space.h"
#include "time_stepping_scheme/00_time_stepping_scheme.h"
#include "time_stepping_scheme/02_time_stepping_scheme_ode.h"
#include "data_management/time_stepping/time_stepping.h"
#include "control/python_config/python_config.h"

namespace ModelOrderReduction
{
  template<typename TimeSteppingType>
  TimeSteppingSchemeOdeReduced<TimeSteppingType>::
  TimeSteppingSchemeOdeReduced(DihuContext context, std::string name):
  MORBase<typename TimeSteppingType::FunctionSpace>(context["ModelOrderReduction"]),
  ::TimeSteppingScheme::TimeSteppingSchemeOdeBase<::FunctionSpace::Generic,1>(context["ModelOrderReduction"],name),
    fullTimestepping_(context["ModelOrderReduction"]), initialized_(false)
  {  
    LOG(DEBUG) << "Constructor TimeSteppingSchemeOdeReduced, given context: " << context.getPythonConfig();

    this->specificSettingsMOR_ = context["ModelOrderReduction"].getPythonConfig();
    
    LOG(DEBUG) << this->specificSettingsMOR_;
    if (this->specificSettingsMOR_.hasKey("nReducedBases"))
    {
      this->nReducedBases_ = this->specificSettingsMOR_.getOptionInt("nReducedBases", 10, PythonUtility::Positive);
      LOG(DEBUG) << "nReducedBases: " << this->nReducedBases_;
    }
    
    if (this->specificSettingsMOR_.hasKey("nRowsSnapshots"))
    {
      this->nRowsSnapshots_ = this->specificSettingsMOR_.getOptionInt("nRowsSnapshots", 10, PythonUtility::Positive);
      LOG(DEBUG) << "nRowsSnapshots: " << this->nRowsSnapshots_;
    }
    
    typedef ::FunctionSpace::Generic GenericFunctionSpace;

    if(this->context_.meshManager()->hasFunctionSpace("functionSpaceReduced"))
    {
      // take the existing function space
      this->functionSpaceRed = this->context_.meshManager()->template functionSpace<GenericFunctionSpace>("functionSpaceReduced");
    }
    else
    {
      // create the functionspace for the reduced order
      LOG(DEBUG) << "nReducedBases: " << this->nReducedBases_;
      this->functionSpaceRed = this->context_.meshManager()->createGenericFunctionSpace(this->nReducedBases_, "functionSpaceReduced");
      LOG(DEBUG) << "functionSpaceRed";
    }
    
    if(this->context_.meshManager()->hasFunctionSpace("functionSpaceRowsSnapshots"))
    {
      // take the existing function space
      this->functionSpaceRowsSnapshots = this->context_.meshManager()->template functionSpace<GenericFunctionSpace>("functionSpaceRowsSnapshots");
    }
    else
    {
      this->functionSpaceRowsSnapshots = this->context_.meshManager()->createGenericFunctionSpace(this->nRowsSnapshots_, "functionSpaceRowsSnapshots");
      LOG(DEBUG) << "functionSpaceRowsSnapshots";
    }
    
    this->data_ = std::make_shared <::Data::TimeStepping<::FunctionSpace::Generic,1>>(context); // create data object
    
  }

  template<typename TimeSteppingType>
  void TimeSteppingSchemeOdeReduced<TimeSteppingType>::setInitialValues()
  {  
    LOG(TRACE) << "TimeSteppingSchemeOdeReduced::setInitialValues()";
    
    Vec &solution = this->fullTimestepping_.data().solution()->getValuesContiguous();
    Vec &redSolution= this->data().solution()->valuesGlobal();
    LOG(DEBUG) << "data solution: " << *this->data().solution();
    
    Mat &basisTransp = this->dataMOR_->basisTransp()->valuesGlobal();
    
    PetscInt mat_sz_1, mat_sz_2;
    PetscInt solution_sz, redSolution_sz;
    
    VecGetSize(solution,&solution_sz);
    VecGetSize(redSolution,&redSolution_sz);
    MatGetSize(basisTransp,&mat_sz_1,&mat_sz_2);
    
    VLOG(2) << "setInitialValues() solution_size: " << solution_sz;
    VLOG(2) << "setInitialValues() redSolution_size: " << redSolution_sz; 
    VLOG(2) << "setInitialValues() mat_sz_1: " << mat_sz_1 << " mat_sz_2: " << mat_sz_2 ;
    
    // reduction step
    this->MatMultReduced(basisTransp, solution, redSolution);

    VLOG(2) << "computed reduced solution: " << *this->data().solution();
  }

  template<typename TimeSteppingType>
  void TimeSteppingSchemeOdeReduced<TimeSteppingType>::
  initialize()
  {
    if (initialized_)
      return;
    
    LOG(TRACE) << "TimeSteppingSchemeOdeReduced::initialize()";

    // add this solver to the solvers diagram, which is an ASCII art representation that will be created at the end of the simulation.
    DihuContext::solverStructureVisualizer()->addSolver("ModelOrderReduction");

    // indicate in solverStructureVisualizer that now a child solver will be initialized
    DihuContext::solverStructureVisualizer()->beginChild("fullTimestepping");

    this->fullTimestepping_.initialize();
    LOG(DEBUG) << "fullTimestepping_ was initialized, has function space: " << this->fullTimestepping_.data().functionSpace()->meshName();

    // indicate in solverStructureVisualizer that the child solver initialization is done
    DihuContext::solverStructureVisualizer()->endChild();

    ::TimeSteppingScheme::TimeSteppingSchemeOdeBase<::FunctionSpace::Generic,1>::initialize();

    this->dataMOR_->setFunctionSpace(this->functionSpaceRed);
    this->dataMOR_->setFunctionSpaceRows(this->functionSpaceRowsSnapshots);
    
    assert(functionSpaceRed->meshPartition());   // assert that the function space was retrieved correctly
    this->data_->setFunctionSpace(functionSpaceRed);
    this->data_->initialize();
    
    MORBase<typename TimeSteppingType::FunctionSpace>::initialize();  
    
    setInitialValues(); //necessary for the explicit scheme
    this->outputWriterManager_.writeOutput(*this->data_, 0, 0);

    VLOG(1) << "initialized full-order solution: " << *this->fullTimestepping_.data().solution();

    LOG(DEBUG) << "fullTimestepping_ has function space: " << this->fullTimestepping_.data().functionSpace()->meshName();

    outputConnectorData_ = std::make_shared<OutputConnectorDataType>();
    outputConnectorData_->addFieldVariable(this->data_->solution());

    // set the outputConnectorData for the solverStructureVisualizer to appear in the solver diagram
    DihuContext::solverStructureVisualizer()->setOutputConnectorData(getOutputConnectorData());

    initialized_ = true;
  }

  template<typename TimeSteppingType>
  void TimeSteppingSchemeOdeReduced<TimeSteppingType>::
  run()
  {
    // initialize
    this->initialize();
    
    // do simulations
    this->advanceTimeSpan();
  }

  template<typename TimeSteppingType>
  TimeSteppingType TimeSteppingSchemeOdeReduced<TimeSteppingType>::
  fullTimestepping()
  {
    return fullTimestepping_;
  }

  template<typename TimeSteppingType>
  std::shared_ptr<typename TimeSteppingSchemeOdeReduced<TimeSteppingType>::OutputConnectorDataType>
  TimeSteppingSchemeOdeReduced<TimeSteppingType>::
  getOutputConnectorData()
  {
    return outputConnectorData_;
  }
} //namespace
