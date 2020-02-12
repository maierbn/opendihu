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

  struct ConfigWithGlobalElementNo
  {
    PythonConfig pythonConfig;      //< python config of the elementWithFaces, containing all settings
    global_no_t elementNoGlobal;    //< global element no, as parsed from config
    element_no_t elementNoLocal;    //< local element no, to be determined
    Mesh::face_t face;              //< face on which the Neumann BC is applied, this is needed here for sorting only
  };

  std::vector<ConfigWithGlobalElementNo> configsWithGlobalElementNo;  //< auxiliary data structure that is used when inputMeshIsGlobal

  // if boundary conditions are present in config
  if (!specificSettings.hasKey(boundaryConditionsConfigKey))
  {
    LOG(WARNING) << specificSettings << "[\"" << boundaryConditionsConfigKey << "\"] not set: No Neumann boundary conditions specified" << std::endl
      << "(Set \"neumannBoundaryConditions\": [] to avoid this warning.)";
    return;
  }

  if (nComponents > 1)
  {
    divideNeumannBoundaryConditionValuesByTotalArea_ = specificSettings.getOptionBool("divideNeumannBoundaryConditionValuesByTotalArea", false);
  }

  // read out if the element nos are specified as local or global values
  bool inputMeshIsGlobal = specificSettings.getOptionBool("inputMeshIsGlobal", true);

  // loop over items in config list
  PyObject *listItem = specificSettings.getOptionListBegin<PyObject*>(boundaryConditionsConfigKey);
  for (;
        !specificSettings.getOptionListEnd(boundaryConditionsConfigKey);
        specificSettings.getOptionListNext<PyObject*>(boundaryConditionsConfigKey, listItem))
  {
    PythonConfig pythonConfigItem = PythonConfig(specificSettings, boundaryConditionsConfigKey, listItem);

    // if element nos are specified as global nos, store them in configsWithGlobalElementNo and transfer to local element nos later
    if (inputMeshIsGlobal)
    {
      // parse necessary information from config object for the ConfigWithGlobalElementNo item
      ConfigWithGlobalElementNo item;
      item.pythonConfig = pythonConfigItem;
      long int elementNoGlobal = pythonConfigItem.getOptionInt("element", 0);

      // negative element nos count from the end
      if (elementNoGlobal < 0)
        elementNoGlobal += functionSpace->nElementsGlobal();

      if (elementNoGlobal < 0 || elementNoGlobal >= functionSpace->nElementsGlobal())
      {
        LOG(ERROR) << "In Neumann boundary conditions, global element no. " << elementNoGlobal << " is invalid (global number of elements: " << functionSpace->nElementsGlobal() << ")";
        continue;
      }
      item.elementNoGlobal = elementNoGlobal;

      std::string faceStr = pythonConfigItem.getOptionString("face", "0+");
      item.face = Mesh::parseFace(faceStr);

      configsWithGlobalElementNo.push_back(item);
    }
    else
    {
      // if element nos are already specified as local nos, directly parse the ElementWithFaces objects
      // parse new Neumann BC, either traction (vector) or flux (scalar), store in boundaryConditionElements_
      boundaryConditionElements_.push_back(parseElementWithFaces(pythonConfigItem, functionSpace, -1));
    }
  }

  // if the element nos were given as global nos, transform them to local nos and parse the ElementWithFaces objects afterwards
  if (inputMeshIsGlobal)
  {
    // sort all configsWithGlobalElementNo by global element no
    auto compareFunction = [](const ConfigWithGlobalElementNo &item1, const ConfigWithGlobalElementNo &item2)
    {
      return item1.elementNoGlobal < item2.elementNoGlobal || (item1.elementNoGlobal == item2.elementNoGlobal && item1.face < item2.face);
    };
    std::sort(configsWithGlobalElementNo.begin(), configsWithGlobalElementNo.end(), compareFunction);

#ifndef NDEBUG
    LOG(DEBUG) << "sorted: ";
    for(ConfigWithGlobalElementNo configWithGlobalElementNo : configsWithGlobalElementNo)
    {
      LOG(DEBUG) << "  elementNoGlobal " << configWithGlobalElementNo.elementNoGlobal << ", face " << Mesh::getString(configWithGlobalElementNo.face);
    }
#endif

    // create copy of configsWithGlobalElementNo with only the locally contained elements in boundaryConditionLocalElements
    std::vector<ConfigWithGlobalElementNo> localConfigsWithGlobalElementNo;

    // loop over all local elements
    for (element_no_t elementNoLocal = 0; elementNoLocal < this->functionSpace_->nElementsLocal(); elementNoLocal++)
    {
      VLOG(2) << "elementNoLocal: " << elementNoLocal;

      // get global element no of the current element
      global_no_t elementNoGlobalNatural = functionSpace_->meshPartition()->getElementNoGlobalNatural(elementNoLocal);

      // check if a boundary condition with this global element no is present in configWithGlobalElementNo
      for (typename std::vector<ConfigWithGlobalElementNo>::const_iterator configBCElementIter = configsWithGlobalElementNo.begin(); configBCElementIter != configsWithGlobalElementNo.end(); configBCElementIter++)
      {
        element_no_t configBCElementNoGlobal = configBCElementIter->elementNoGlobal;   // this is now the stored global no, because inputMeshIsGlobal

        if (configBCElementNoGlobal == elementNoGlobalNatural)
        {
          // add to boundaryConditionLocalElements
          localConfigsWithGlobalElementNo.push_back(*configBCElementIter);

          // set element no to real local no
          localConfigsWithGlobalElementNo.back().elementNoLocal = elementNoLocal;
          break;
        }
        else if (configBCElementNoGlobal > elementNoGlobalNatural)
        {
          break;
        }
      }
    }

    // assign contents of boundaryConditionLocalElements back to boundaryConditionElements_
    for (const ConfigWithGlobalElementNo &configWithGlobalElementNo : localConfigsWithGlobalElementNo)
    {
      // parse new Neumann BC, either traction (vector) or flux (scalar), store in boundaryConditionElements_
      boundaryConditionElements_.push_back(parseElementWithFaces(configWithGlobalElementNo.pythonConfig, functionSpace, configWithGlobalElementNo.elementNoLocal));
    }
  }

  // output parsed settings
  LOG(DEBUG) << boundaryConditionElements_.size() << " boundaryConditionElements parsed";
  LOG(DEBUG) << boundaryConditionElements_.size() << " boundaryConditionElements_: ";

  for(ElementWithFaces elementWithFaces : boundaryConditionElements_)
  {
    LOG(DEBUG) << "  elementNoLocal " << elementWithFaces.elementNoLocal << ", face " << Mesh::getString(elementWithFaces.face)
      << ", dofVectors on surface: " << elementWithFaces.dofVectors << ", surfaceDofs on volume: " << elementWithFaces.surfaceDofs;
  }

  initializeRhs();
}

// initialize directly
template<typename FunctionSpaceType,typename QuadratureType,int nComponents>
void NeumannBoundaryConditionsBase<FunctionSpaceType,QuadratureType,nComponents>::
initialize(std::shared_ptr<FunctionSpaceType> functionSpace, const std::vector<typename NeumannBoundaryConditionsBase<FunctionSpaceType,QuadratureType,nComponents>::ElementWithFaces> &boundaryConditionElements)
{
  this->functionSpace_ = functionSpace;
  data_.setFunctionSpace(functionSpace_);
  data_.initialize();

  boundaryConditionElements_ = boundaryConditionElements;

  // output parsed settings
  LOG(DEBUG) << boundaryConditionElements_.size() << " boundaryConditionElements parsed";
  LOG(DEBUG) << boundaryConditionElements_.size() << " boundaryConditionElements_: ";

  for(ElementWithFaces elementWithFaces : boundaryConditionElements_)
  {
    LOG(DEBUG) << "  elementNoLocal " << elementWithFaces.elementNoLocal << ", face " << Mesh::getString(elementWithFaces.face)
      << ", dofVectors on surface: " << elementWithFaces.dofVectors << ", surfaceDofs on volume: " << elementWithFaces.surfaceDofs;
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
