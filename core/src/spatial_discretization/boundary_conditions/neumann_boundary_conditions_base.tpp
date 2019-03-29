#include "spatial_discretization/boundary_conditions/neumann_boundary_conditions_base.h"

#include "easylogging++.h"
#include "utility/python_utility.h"
#include "utility/vector_operators.h"
#include "control/types.h"

namespace SpatialDiscretization
{

//constructor
template<typename FunctionSpaceType,typename QuadratureType,int nComponents>
NeumannBoundaryConditionsBase<FunctionSpaceType,QuadratureType,nComponents>::NeumannBoundaryConditionsBase(DihuContext context) :
  BoundaryConditionsBase<FunctionSpaceType, nComponents>::BoundaryConditionsBase(context), data_(context)
{
}

// set the boundary conditions to the right hand side
template<typename FunctionSpaceType,typename QuadratureType,int nComponents>
void NeumannBoundaryConditionsBase<FunctionSpaceType,QuadratureType,nComponents>::
initialize(PythonConfig specificSettings, std::shared_ptr<FunctionSpaceType> functionSpace, std::string boundaryConditionsConfigKey)
{
  this->functionSpace_ = functionSpace;
  data_.setFunctionSpace(functionSpace_);
  data_.initialize();

  LOG(DEBUG) << "initialize Neumann boundary conditions, key  \"" << boundaryConditionsConfigKey << "\"";

  // if boundary conditions are present in config
  if (specificSettings.hasKey(boundaryConditionsConfigKey))
  {
    // loop over items in config list
    PyObject *listItem = specificSettings.getOptionListBegin<PyObject*>(boundaryConditionsConfigKey);
    for (;
         !specificSettings.getOptionListEnd(boundaryConditionsConfigKey);
         specificSettings.getOptionListNext<PyObject*>(boundaryConditionsConfigKey, listItem))
    {
      // parse new Neumann BC, either traction (vector) or flux (scalar), store in boundaryConditionElements_
      boundaryConditionElements_.push_back(parseElementWithFaces(listItem, functionSpace));
    }
  }
  else
  {
    LOG(WARNING) << specificSettings << "[\"" << boundaryConditionsConfigKey << "\"] not set: No Neumann boundary conditions specified";
  }

  LOG(DEBUG) << boundaryConditionElements_.size() << " boundaryConditionElements parsed";

  // if element nos are given as global nos, set to local element nos

  // determine if the BC indices in the config are given for global or local dof nos
  bool inputMeshIsGlobal = specificSettings.getOptionBool("inputMeshIsGlobal", true);
  if (inputMeshIsGlobal)
  {
    // sort all boundaryConditionElements_ by (here global) element no
    auto compareFunction = [](const ElementWithFaces &item1, const ElementWithFaces &item2)
    {
      return item1.elementNoLocal < item2.elementNoLocal || (item1.elementNoLocal == item2.elementNoLocal && item1.face < item2.face);
    };
    std::sort(boundaryConditionElements_.begin(), boundaryConditionElements_.end(), compareFunction);

    LOG(DEBUG) << "sorted: ";

    for(ElementWithFaces elementWithFaces : boundaryConditionElements_)
    {
      LOG(DEBUG) << "  elementNoLocal " << elementWithFaces.elementNoLocal << ", face " << Mesh::getString(elementWithFaces.face) << ", dofVectors: " << elementWithFaces.dofVectors;
    }

    // create copy of boundaryConditionElements_ with only the locally contained elements in boundaryConditionLocalElements
    std::vector<ElementWithFaces> boundaryConditionLocalElements;

    // loop over all local elements
    for (element_no_t elementNoLocal = 0; elementNoLocal < this->functionSpace_->nElementsLocal(); elementNoLocal++)
    {
      VLOG(2) << "elementNoLocal: " << elementNoLocal;

      // get global element no
      global_no_t elementNoGlobalNatural = functionSpace_->meshPartition()->getElementNoGlobalNatural(elementNoLocal);

      // check if a boundary condition with this global element no is present
      for (typename std::vector<ElementWithFaces>::const_iterator configBCElementIter = boundaryConditionElements_.begin(); configBCElementIter != boundaryConditionElements_.end(); configBCElementIter++)
      {
        element_no_t configBCElementNoGlobal = configBCElementIter->elementNoLocal;   // this is now the global no, because inputMeshIsGlobal

        if (configBCElementNoGlobal == elementNoGlobalNatural)
        {
          // add to boundaryConditionLocalElements
          boundaryConditionLocalElements.push_back(*configBCElementIter);
          // set element no to real local no
          boundaryConditionLocalElements.back().elementNoLocal = elementNoLocal;
          break;
        }
        else if (configBCElementNoGlobal > elementNoGlobalNatural)
        {
          break;
        }
      }
    }

    // assign contents of boundaryConditionLocalElements back to boundaryConditionElements_
    boundaryConditionElements_.assign(boundaryConditionLocalElements.begin(), boundaryConditionLocalElements.end());
  }

  LOG(DEBUG) << boundaryConditionElements_.size() << " boundaryConditionElements_: ";

  for(ElementWithFaces elementWithFaces : boundaryConditionElements_)
  {
    LOG(DEBUG) << "  elementNoLocal " << elementWithFaces.elementNoLocal << ", face " << Mesh::getString(elementWithFaces.face) << ", dofVectors: " << elementWithFaces.dofVectors;
  }

  initializeRhs();
}

template<typename FunctionSpaceType,typename QuadratureType,int nComponents>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> NeumannBoundaryConditionsBase<FunctionSpaceType,QuadratureType,nComponents>::
rhs()
{
  return data_.rhs();
}

}  // namespace
