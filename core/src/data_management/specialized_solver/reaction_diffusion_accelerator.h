#pragma once

#include <Python.h>  // has to be the first included header
/*#include <petscmat.h>
#include <memory>

#include "control/types.h"
#include "mesh/mesh.h"
#include "data_management/data.h"
#include "data_management/time_stepping/time_stepping.h"
#include "field_variable/field_variable.h"*/

namespace Data
{

/**  The datastructures used for static bidomain solver.
  */
template<typename FunctionSpaceType,int nStates>
class ReactionDiffusionAccelerator : public Data<FunctionSpaceType>
{
public:

  typedef FieldVariable::FieldVariable<FunctionSpaceType,nStates> StatesFieldVariableType;
  typedef FieldVariable::FieldVariable<FunctionSpaceType,1> IntermediatesFieldVariableType;
  typedef FieldVariable::FieldVariable<FunctionSpaceType,2> ParametersFieldVariableType;
  typedef FieldVariable::FieldVariable<FunctionSpaceType,3> GeometryFieldVariableType;

  //! constructor
  ReactionDiffusionAccelerator(DihuContext context);

  //! return a reference to the states vector
  std::shared_ptr<StatesFieldVariableType> states();

  //! return the field variable of the rates
  std::shared_ptr<StatesFieldVariableType> rates();

  //! return a reference to the intermediates vector
  std::shared_ptr<IntermediatesFieldVariableType> intermediates();

  //! return a reference to the parameters vector
  std::shared_ptr<ParametersFieldVariableType> parameters();

  //! initialize
  void initialize();

  //! field variables that will be output by outputWriters
  typedef std::tuple<
    std::shared_ptr<GeometryFieldVariableType>,     // geometry
    std::shared_ptr<StatesFieldVariableType>,     // states
    std::shared_ptr<StatesFieldVariableType>,     // rates
    std::shared_ptr<IntermediatesFieldVariableType>,              // intermediates
    std::shared_ptr<ParametersFieldVariableType>              // parameters
  > OutputFieldVariables;

  //! get pointers to all field variables that can be written by output writers
  OutputFieldVariables getOutputFieldVariables();

private:

  //! initializes the vectors with size
  void createPetscObjects() override;

  std::shared_ptr<StatesFieldVariableType> states_; ///< the states of all cells, usually 57 components each (for Shorten)
  std::shared_ptr<StatesFieldVariableType> rates_; ///< the rates (d states / dt ) of all cells, usually 57 components each
  std::shared_ptr<IntermediatesFieldVariableType> intermediates_; ///< values of interes needed outsid of 0D
  std::shared_ptr<ParametersFieldVariableType> parameters_; ///<  values entering 0D as parameter. 
  // std::shared_ptr<FieldVariableType> constants_; ///< 
  // std::shared_ptr<FieldVariableType> vecToActualizeDiffusionCells_;
  // std::shared_ptr<FieldVariableType> interimResultsVector_;
};

} // namespace Data

#include "data_management/specialized_solver/reaction_diffusion_accelerator.tpp"
