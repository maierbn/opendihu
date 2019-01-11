#include "operator_splitting/solution_vector_mapping/solution_vector_mapping.h"

#include <vector>
#include <tuple>
#include "easylogging++.h"

template<typename BasisFunctionType, typename FieldVariableType2>
void SolutionVectorMapping<
  std::vector<std::vector<std::tuple<std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType>>, int, double>>>,   // vector<vector<fieldVariableType,componentNo,prefactor>>
  std::shared_ptr<FieldVariableType2>  // <3D field variable>
>::transfer(const std::vector<std::vector<std::tuple<std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType>>, int, double>>> &transferableSolutionData1,
            std::shared_ptr<FieldVariableType2> transferableSolutionData2)
{
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType> FieldVariableType1;

  // loop over vector of fibers
  for (int i = 0; i < transferableSolutionData1.size(); i++)
  {
    for (int j = 0; j < transferableSolutionData1[i].size(); j++)
    {
       std::shared_ptr<FieldVariableType1> transmembranePotential = std::get<0>(transferableSolutionData1[i][j]);
       std::shared_ptr<MappingBetweenMeshes<FieldVariableType1::FunctionSpace, FieldVariableType2::FunctionSpace>> mappingBetweenMeshes
         = DihuContext::meshManager()->mappingBetweenMeshes(transmembranePotential->functionSpace()->meshName(), transferableSolutionData2->functionSpace()->meshName());
       mappingBetweenMeshes->map(transmembranePotential, 0, transferableSolutionData2, 0);
    }
  }
};

template<typename BasisFunctionType, typename FieldVariableType1>
void SolutionVectorMapping<
  std::shared_ptr<FieldVariableType1>,  // <3D field variable>
  std::vector<std::vector<std::tuple<std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType>>, int, double>>>   // vector<vector<fieldVariableType,componentNo,prefactor>>
>::transfer(std::shared_ptr<FieldVariableType2> transferableSolutionData1,
            std::vector<std::vector<std::tuple<std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunctionType>>, int, double>>> &transferableSolutionData2)
{

};

template<typename FieldVariableType1, typename FieldVariableType2>
void SolutionVectorMapping<
  std::vector<std::tuple<std::shared_ptr<FieldVariableType1>, int, double>>,   // <fieldVariableType,componentNo,prefactor>
  std::pair<std::vector<Vec>,std::vector<std::shared_ptr<FieldVariableType2>>>  // Petsc Vecs which are the sub-vectors of a nested vector
>::transfer(const std::vector<std::tuple<std::shared_ptr<FieldVariableType1>, int, double>> &transferableSolutionData1,
            std::pair<std::vector<Vec>,std::vector<std::shared_ptr<FieldVariableType2>>> transferableSolutionData2)
{
  // transfer from cellml to multidomain
  // this is analogous to the following code line in the old code:
  // dataMultidomain_.subcellularStates(k)->extractComponentCopy(0, dataMultidomain_.transmembranePotential(k));

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

template<typename FieldVariableType1, typename FieldVariableType2>
void SolutionVectorMapping<
  std::pair<std::vector<Vec>,std::vector<std::shared_ptr<FieldVariableType1>>>,  // Petsc Vecs which are the sub-vectors of a nested vector
  std::vector<std::tuple<std::shared_ptr<FieldVariableType2>, int, double>>   // <fieldVariable,componentNo,prefactor>
>::transfer(const std::pair<std::vector<Vec>,std::vector<std::shared_ptr<FieldVariableType1>>> &transferableSolutionData1,  // <Petsc Vecs which are the sub-vectors of a nested vector, transmembranePotential>
            const std::vector<std::tuple<std::shared_ptr<FieldVariableType2>, int, double>> &transferableSolutionData2)
{
  // transfer from multidomain to cellml
  // this is analogous to the following code line in the old code:
  // dataMultidomain_.subcellularStates(k)->setValues(0, subVectors[k]);   // note, subcellularStates is in contiguous representation

  // rename input variables
  const std::vector<Vec> &subVectors = transferableSolutionData1.first;
  const std::vector<std::shared_ptr<FieldVariableType1>> &transmembranePotential = transferableSolutionData1.second;
  const std::vector<std::tuple<std::shared_ptr<FieldVariableType2>, int, double>> &subcellularStates = transferableSolutionData2;

  // assert that dimensions match
  int nCompartments = subVectors.size() - 1;  // not for phi_e
  LOG(DEBUG) << "transfer from multidomain to cellml, size of subvectors=" << nCompartments << ", subcellularStates.size(): " << subcellularStates.size();
  assert(nCompartments == subcellularStates.size());

  // transfer from subVectors[k] to subcellularStatesFieldVariable
  // loop over compartments and transfer data from subVectors of multidomain solver to subcellularState(0) (=Vm)
  for (int k = 0; k < nCompartments; k++)
  {
    std::shared_ptr<FieldVariableType2> subcellularStatesFieldVariable = std::get<0>(subcellularStates[k]);
    int outputStateIndex = std::get<1>(subcellularStates[k]);  // is 0

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

/*
/store/software/opendihu/core/src/operator_splitting/coupling_or_godunov.tpp:53:15: error: ‘transfer’ is not a member of ‘
SolutionVectorMapping<
  std::vector<
    std::vector<
      std::tuple<
        std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1> >, 1> >,
        int,
        double
      >,
      std::allocator<
        std::tuple<
          std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1> >, 1> >,
          int, double>
      >
    >,
    std::allocator<
      std::vector<
        std::tuple<
          std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1> >, 1> >, int, double
        >,
        std::allocator<
          std::tuple<
            std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1> >, 1> >, int, double
          >
        >
      >
    >
  >,

  std::pair<
    std::vector<_p_Vec*>,
    std::vector<
      std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunction::LagrangeOfOrder<1> >, 1> >, std::allocator<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunction::LagrangeOfOrder<1> >, 1> > >
    >
  >
>
’
 */
