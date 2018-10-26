#include "operator_splitting/solution_vector_mapping/solution_vector_mapping.h"

#include <vector>
#include <tuple>
#include "easylogging++.h"

template<typename FieldVariableType1, typename FieldVariableType2>
void SolutionVectorMapping<
  std::vector<std::tuple<std::shared_ptr<FieldVariableType1>, int, double>>,   // <fieldVariableType,componentNo,prefactor>
  std::pair<std::vector<Vec>,std::vector<std::shared_ptr<FieldVariableType2>>>  // Petsc Vecs which are the sub-vectors of a nested vector
>::transfer(const std::vector<std::tuple<std::shared_ptr<FieldVariableType1>, int, double>> &transferableSolutionData1,
            std::pair<std::vector<Vec>,std::vector<std::shared_ptr<FieldVariableType2>>> transferableSolutionData2)
{
  // transfer from cellml to multidomain
  // this is analogous to the following code line in the old code:
  // dataMultidomain_.subcellularStates(k)->extractComponent(0, dataMultidomain_.transmembranePotential(k));

  // rename input variables
  const std::vector<std::tuple<std::shared_ptr<FieldVariableType1>, int, double>> &subcellularStates = transferableSolutionData1;
  std::vector<std::shared_ptr<FieldVariableType2>> &transmembranePotential = transferableSolutionData2.second;

  // assert that dimensions match
  int nCompartments = transmembranePotential.size();
  LOG(DEBUG) << "transfer from cellml to multidomain, nCompartments=" << nCompartments << ", subcellularStates.size(): " << subcellularStates.size();
  assert(nCompartments == subcellularStates.size());

  // loop over compartments and transfer data from subcellularState(0) (=Vm) to transmembranePotential of multidomain solver
  for (int k = 0; k < nCompartments; k++)
  {
    std::shared_ptr<FieldVariableType1> subcellularStatesFieldVariable = std::get<0>(subcellularStates[k]);
    int outputStateIndex = std::get<1>(subcellularStates[k]);  // is 0
    double prefactor = std::get<2>(subcellularStates[k]);  // is 1
    std::shared_ptr<FieldVariableType2> transmembranePotentialCompartment = transmembranePotential[k];
    subcellularStatesFieldVariable->extractComponent(outputStateIndex, transmembranePotentialCompartment);

    if (prefactor != 1.0)
    {
      PetscErrorCode ierr;
      ierr = VecScale(transmembranePotential[k]->valuesGlobal(), prefactor); CHKERRV(ierr);
    }

  }
}

template<typename FieldVariableType1, typename FieldVariableType2>
void SolutionVectorMapping<
  std::pair<std::vector<Vec>,std::vector<std::shared_ptr<FieldVariableType1>>>,  // Petsc Vecs which are the sub-vectors of a nested vector
  std::vector<std::tuple<std::shared_ptr<FieldVariableType2>, int, double>>   // <fieldVariable,componentNo,prefactor>
>::transfer(const std::pair<std::vector<Vec>,std::vector<std::shared_ptr<FieldVariableType1>>> &transferableSolutionData1,
            const std::vector<std::tuple<std::shared_ptr<FieldVariableType2>, int, double>> &transferableSolutionData2)
{
  // transfer from multidomain to cellml
  // this is analogous to the following code line in the old code:
  // dataMultidomain_.subcellularStates(k)->setValues(0, subVectors[k]);   // note, subcellularStates is in contiguous representation

  // rename input variables
  const std::vector<Vec> &subVectors = transferableSolutionData1.first;
  const std::vector<std::tuple<std::shared_ptr<FieldVariableType2>, int, double>> &subcellularStates = transferableSolutionData2;

  // assert that dimensions match
  int nCompartments = subVectors.size() - 1;  // not for phi_e
  LOG(DEBUG) << "transfer from multidomain to cellml, size of subvectors=" << nCompartments << ", subcellularStates.size(): " << subcellularStates.size();
  assert(nCompartments == subcellularStates.size());

  // loop over compartments and transfer data from subVectors of multidomain solver to subcellularState(0) (=Vm)
  for (int k = 0; k < nCompartments; k++)
  {

    std::shared_ptr<FieldVariableType2> subcellularStatesFieldVariable = std::get<0>(subcellularStates[k]);
    int outputStateIndex = std::get<1>(subcellularStates[k]);  // is 0

    subcellularStatesFieldVariable->setValues(outputStateIndex, subVectors[k]);
  }
}
