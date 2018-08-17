#include "output_writer/python/python.h"

#include <Python.h>   // this has to be the first included header
#include <iostream>

#include "easylogging++.h"
#include "utility/python_utility.h"

namespace OutputWriter
{

template<int D, typename BasisFunctionType, typename OutputFieldVariablesType>
PyObject *Python<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,OutputFieldVariablesType>::
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
  typedef BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType> BasisOnMeshType;
  std::shared_ptr<BasisOnMeshType> mesh = std::static_pointer_cast<BasisOnMeshType>(meshBase);

  std::string basisFunction = BasisOnMeshType::BasisFunction::getBasisFunctionString();
  int basisOrder = BasisOnMeshType::BasisFunction::getBasisOrder();

  // start critical section for python API calls
  PythonUtility::GlobalInterpreterLock lock;
  
  PyObject *pyElementalDofs = Python<BasisOnMeshType,OutputFieldVariablesType>::
    buildPyElementalDofsObject(meshBase, onlyNodalValues);
  
  // build python dict that will contain all information and data
  PyObject *data = Py_BuildValue("{s s, s i, s i, s s, s i, s O, s O, s O, s i, s d}", "meshType", "UnstructuredDeformable",
                                 "dimension", D, "nElements", mesh->nElementsLocal(),
                                 "basisFunction", basisFunction.c_str(), "basisOrder", basisOrder,
                                 "onlyNodalValues", onlyNodalValues ? Py_True: Py_False,
                                 "data", pyData, "elementalDofs", pyElementalDofs, 
                                 "timeStepNo", timeStepNo, "currentTime", currentTime);

  //LOG(DEBUG) << data << " done";

  return data;
}

template<int D, typename BasisFunctionType, typename OutputFieldVariablesType>
PyObject *Python<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,OutputFieldVariablesType>::
buildPyElementalDofsObject(std::shared_ptr<Mesh::Mesh> meshBase, bool onlyNodalValues)
{
  typedef BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType> BasisOnMesh;
  std::shared_ptr<BasisOnMesh> mesh = std::static_pointer_cast<BasisOnMesh>(meshBase);

  // create a list of lists for each element the node numbers (if onlyNodalValues) or the dofs
  PyObject *pyElementalDofs = PyList_New((Py_ssize_t)mesh->nElementsLocal());
  
  // loop over elements
  for (element_no_t elementNo = 0; elementNo < mesh->nElementsLocal(); elementNo++)
  {
    std::vector<node_no_t> dofs;
    
    std::array<dof_no_t,BasisOnMesh::nDofsPerElement()> dofsOfElement = mesh->getElementDofLocalNos(elementNo);
    for (typename std::array<dof_no_t,BasisOnMesh::nDofsPerElement()>::const_iterator iter = dofsOfElement.begin(); iter != dofsOfElement.end(); iter++)
    {
      dof_no_t dofNo = *iter;
      
      if (onlyNodalValues)
      {
        if (dofNo % BasisOnMesh::nDofsPerNode() == 0)
        {
          node_no_t nodeNo = dofNo / BasisOnMesh::nDofsPerNode();
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
};