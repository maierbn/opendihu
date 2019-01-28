#include "output_writer/python/python.h"

#include <Python.h>   // this has to be the first included header
#include <iostream>

#include "easylogging++.h"
#include "utility/python_utility.h"

namespace OutputWriter
{

template<int D, typename BasisFunctionType, typename OutputFieldVariablesType>
PyObject *Python<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,OutputFieldVariablesType>::
buildPyDataObject(OutputFieldVariablesType fieldVariables, 
                  std::string meshName, int timeStepNo, double currentTime, bool onlyNodalValues)
{
  // build python dict containing all information
  // data = {
  //   "meshType" : "UnstructuredDeformable",
  //   "dimension": dim,
  //   "nElements" : nElements,
  //   "basisFunction" : "Lagrange",
  //   "basisOrder" : "1",
  //   "onlyNodalValues" : True,
  //   "data" : [
  //      {"name" : "fieldVariableName",
  //       "components" : [
  //           {"name" : "componentName", "values": data},
  //       ]
  //      },
  //   ]
  //   "timeStepNo" : timeStepNo,
  //   "currentTime" : currentTime
  // }


  LOG(DEBUG) << "build pyData, " << std::tuple_size<OutputFieldVariablesType>::value;

  // build python object for data
  std::shared_ptr<Mesh::Mesh> meshBase;
  PyObject *pyData = PythonBase<OutputFieldVariablesType>::buildPyFieldVariablesObject(fieldVariables, meshName, onlyNodalValues, meshBase);

  // cast mesh to its real type
  typedef FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType> FunctionSpaceType;
  std::shared_ptr<FunctionSpaceType> mesh = std::static_pointer_cast<FunctionSpaceType>(meshBase);

  std::string basisFunction = FunctionSpaceType::BasisFunction::getBasisFunctionString();
  int basisOrder = FunctionSpaceType::BasisFunction::getBasisOrder();

  int nRanks = mesh->meshPartition()->nRanks();
  int ownRankNo = mesh->meshPartition()->ownRankNo();

  // start critical section for python API calls
  // PythonUtility::GlobalInterpreterLock lock;
  
  PyObject *pyElementalDofs = Python<FunctionSpaceType,OutputFieldVariablesType>::
    buildPyElementalDofsObject(meshBase, onlyNodalValues);
  
  // build python dict that will contain all information and data
  PyObject *data = Py_BuildValue("{s s, s s, s s, s i, s i, s s, s i, s O, s i, s i, s O, s O, s i, s d}",
                                 "version", DihuContext::versionText().c_str(),
                                 "meta", DihuContext::metaText().c_str(),
                                 "meshType", "UnstructuredDeformable",
                                 "dimension", D, "nElements", mesh->nElementsLocal(),
                                 "basisFunction", basisFunction.c_str(), "basisOrder", basisOrder,
                                 "onlyNodalValues", onlyNodalValues ? Py_True: Py_False,
                                 "nRanks", nRanks, "ownRankNo", ownRankNo,
                                 "data", pyData, "elementalDofs", pyElementalDofs, 
                                 "timeStepNo", timeStepNo, "currentTime", currentTime);

  //LOG(DEBUG) << data << " done";

  return data;
}

template<int D, typename BasisFunctionType, typename OutputFieldVariablesType>
PyObject *Python<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,OutputFieldVariablesType>::
buildPyElementalDofsObject(std::shared_ptr<Mesh::Mesh> meshBase, bool onlyNodalValues)
{
  typedef FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType> FunctionSpace;
  std::shared_ptr<FunctionSpace> mesh = std::static_pointer_cast<FunctionSpace>(meshBase);

  // create a list of lists for each element the node numbers (if onlyNodalValues) or the dofs
  PyObject *pyElementalDofs = PyList_New((Py_ssize_t)mesh->nElementsLocal());
  
  // loop over elements
  for (element_no_t elementNo = 0; elementNo < mesh->nElementsLocal(); elementNo++)
  {
    std::vector<node_no_t> dofs;
    
    std::array<dof_no_t,FunctionSpace::nDofsPerElement()> dofsOfElement = mesh->getElementDofNosLocal(elementNo);
    for (typename std::array<dof_no_t,FunctionSpace::nDofsPerElement()>::const_iterator iter = dofsOfElement.begin(); iter != dofsOfElement.end(); iter++)
    {
      dof_no_t dofNo = *iter;
      
      if (onlyNodalValues)
      {
        if (dofNo % FunctionSpace::nDofsPerNode() == 0)
        {
          node_no_t nodeNo = dofNo / FunctionSpace::nDofsPerNode();
          dofs.push_back(nodeNo);
        }
      }
      else 
      {
        dofs.push_back(dofNo);
      }
    }
    
    PyObject *pyDofs = PythonUtility::convertToPythonList(dofs);
    
    // add to list
    PyList_SetItem(pyElementalDofs, (Py_ssize_t)elementNo, pyDofs);    // steals reference to pyComponent
  }
  
  return pyElementalDofs;
}
}  // namespace
