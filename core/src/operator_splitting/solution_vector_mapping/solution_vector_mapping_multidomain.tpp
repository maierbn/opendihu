#include "operator_splitting/solution_vector_mapping/solution_vector_mapping.h"

#include <vector>
#include <tuple>
#include "easylogging++.h"

template<typename FunctionSpaceType1, int nComponents1, typename FieldVariableType2>
void SolutionVectorMapping<
  std::vector<Data::ScaledFieldVariableComponent<FunctionSpaceType1,nComponents1>>,   // <fieldVariableType,componentNo,prefactor>
  std::pair<std::vector<Vec>,std::vector<std::shared_ptr<FieldVariableType2>>>  // Petsc Vecs which are the sub-vectors of a nested vector
>::transfer(const std::vector<Data::ScaledFieldVariableComponent<FunctionSpaceType1,nComponents1>> &transferableSolutionData1,
            std::pair<std::vector<Vec>,std::vector<std::shared_ptr<FieldVariableType2>>> transferableSolutionData2,
            const std::string transferSlotName)
{
  // transfer from cellml to multidomain
  // this is analogous to the following code line in the old code:
  // dataMultidomain_.subcellularStates(k)->extractComponentCopy(0, dataMultidomain_.transmembranePotential(k));

  // rename input variables
  const std::vector<Data::ScaledFieldVariableComponent<FunctionSpaceType1,nComponents1>> &subcellularStates = transferableSolutionData1;
  std::vector<std::shared_ptr<FieldVariableType2>> &transmembranePotential = transferableSolutionData2.second;
  typedef FieldVariable::FieldVariable<FunctionSpaceType1,nComponents1> FieldVariableType1;

  // assert that dimensions match
  int nCompartments = transmembranePotential.size();
  VLOG(1) << "Solution vector mapping (solution_vector_mapping_multidomain.tpp)";
  LOG(DEBUG) << "transfer from cellml to multidomain, nCompartments=" << nCompartments << ", subcellularStates.size(): " << subcellularStates.size();
  assert(nCompartments == subcellularStates.size());

  // loop over compartments and transfer data from subcellularState(0) (=Vm) to transmembranePotential of multidomain solver
  for (int k = 0; k < nCompartments; k++)
  {
    std::shared_ptr<FieldVariableType1> subcellularStatesFieldVariable = subcellularStates[k].values;
    int outputStateIndex = subcellularStates[k].componentNo;  // is 0
    double prefactor = subcellularStates[k].scalingFactor;  // is 1
    std::shared_ptr<FieldVariableType2> transmembranePotentialCompartment = transmembranePotential[k];
    //subcellularStatesFieldVariable->extractComponentShared(outputStateIndex, transmembranePotentialCompartment);  // extract memory location, here it is possible to reuse data
    subcellularStatesFieldVariable->extractComponentCopy(outputStateIndex, transmembranePotentialCompartment);      // copy vector

    if (prefactor != 1.0)
    {
      PetscErrorCode ierr;
      ierr = VecScale(transmembranePotential[k]->valuesGlobal(), prefactor); CHKERRV(ierr);
    }
    LOG(DEBUG) << "transmembranePotential[" << k << "]: " << *transmembranePotential[k];
  }
}

template<typename FieldVariableType1, typename FunctionSpaceType2, int nComponents2>
void SolutionVectorMapping<
  std::pair<std::vector<Vec>,std::vector<std::shared_ptr<FieldVariableType1>>>,  // Petsc Vecs which are the sub-vectors of a nested vector
  std::vector<Data::ScaledFieldVariableComponent<FunctionSpaceType2,nComponents2>>   // <fieldVariable,componentNo,prefactor>
>::transfer(const std::pair<std::vector<Vec>,std::vector<std::shared_ptr<FieldVariableType1>>> &transferableSolutionData1,  // <Petsc Vecs which are the sub-vectors of a nested vector, transmembranePotential>
            const std::vector<Data::ScaledFieldVariableComponent<FunctionSpaceType2,nComponents2>> &transferableSolutionData2,
            const std::string transferSlotName)
{
  // transfer from multidomain to cellml
  // this is analogous to the following code line in the old code:
  // dataMultidomain_.subcellularStates(k)->setValues(0, subVectors[k]);   // note, subcellularStates is in contiguous representation

  // rename input variables
  const std::vector<Vec> &subVectors = transferableSolutionData1.first;
  const std::vector<std::shared_ptr<FieldVariableType1>> &transmembranePotential = transferableSolutionData1.second;
  const std::vector<Data::ScaledFieldVariableComponent<FunctionSpaceType2,nComponents2>> &subcellularStates = transferableSolutionData2;
  typedef FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2> FieldVariableType2;

  // assert that dimensions match
  int nCompartments = subVectors.size() - 1;  // not for phi_e
  VLOG(1) << "Solution vector mapping (solution_vector_mapping_multidomain.tpp)";
  LOG(DEBUG) << "transfer from multidomain to cellml, size of subvectors=" << nCompartments << ", subcellularStates.size(): " << subcellularStates.size();
  assert(nCompartments == subcellularStates.size());

  // transfer from subVectors[k] to subcellularStatesFieldVariable
  // loop over compartments and transfer data from subVectors of multidomain solver to subcellularState(0) (=Vm)
  for (int k = 0; k < nCompartments; k++)
  {
    std::shared_ptr<FieldVariableType2> subcellularStatesFieldVariable = subcellularStates[k].values;
    int outputStateIndex = subcellularStates[k].componentNo;  // is 0

    // restore from subcellularStatesFieldVariable extracted data vector to allow usage of subcellularStatesFieldVariable
    if (subcellularStatesFieldVariable->partitionedPetscVec()->currentRepresentation() == Partition::values_representation_t::representationInvalid)
    {
      subcellularStatesFieldVariable->restoreExtractedComponent(transmembranePotential[k]->partitionedPetscVec());
    }

    // copy the values to their location in subcellularStatesFieldVariable, this copy cannot be avoided
    subcellularStatesFieldVariable->setRepresentationGlobal();
    subcellularStatesFieldVariable->setValues(outputStateIndex, subVectors[k]);

    LOG(DEBUG) << "subcellularStates[" << k << "]: " << *subcellularStatesFieldVariable;
  }
}
