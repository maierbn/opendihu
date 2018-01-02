#include "output_writer/callback.h"

#include <iostream>

#include "easylogging++.h"
#include <Python.h>

#include <utility/python_utility.h>
#include <utility/petsc_utility.h>
#include <mesh/regular_fixed.h>
#include <mesh/structured_deformable.h>
#include <mesh/mesh.h>

namespace OutputWriter
{

Callback::Callback(PyObject *settings) : Generic(settings)
{
  callback_ = PythonUtility::getOptionPyObject(settings, "callback");
}

void Callback::callCallback(std::vector<double> &data, std::vector<long> &nEntries, int timeStepNo, double currentTime)
{
  if (callback_ == NULL)
    return;
  
  // compose callback function
  PyObject *dataList = PythonUtility::convertToPythonList(data);
  PyObject *nEntriesList = PythonUtility::convertToPythonList(nEntries);
  
  // signature: def callback(data, shape, nEntries, dimension, timeStepNo, currentTime)
  PyObject *arglist = Py_BuildValue("(O,O,i,i,i,d)", dataList, nEntriesList, data.size(), nEntries.size(), timeStepNo, currentTime);
  PyObject *returnValue = PyObject_CallObject(callback_, arglist);
  
  // if there was an error while executing the function, print the error message
  if (returnValue == NULL)
    PyErr_Print();
  
  // decrement reference counters for python objects
  Py_DECREF(dataList);
  Py_DECREF(nEntriesList);
  Py_DECREF(returnValue);
  Py_DECREF(arglist);
}

void Callback::writeSolution(Data::Data& data, int timeStepNo, double currentTime)
{
  if (!data.mesh())
  {
    LOG(FATAL) << "Mesh is not set!";
  }
  const int dimension = data.mesh()->dimension();
  
  // solution and rhs vectors in mesh shape
  switch(dimension)
  {
  case 1:
    writeSolutionDim<1>(data, timeStepNo, currentTime);
    break;
  case 2:
    writeSolutionDim<2>(data, timeStepNo, currentTime);
    break;
  case 3:
    writeSolutionDim<3>(data, timeStepNo, currentTime);
    break;
  };
}

template <int dimension>
void Callback::writeSolutionDim(Data::Data &data, int timeStepNo, double currentTime)
{
  // solution and rhs vectors in mesh shape
  // if mesh is regular fixed
  if (std::dynamic_pointer_cast<Mesh::RegularFixed<dimension>>(data.mesh()) != NULL)
  {
    // cast mesh to its special type
    std::shared_ptr<Mesh::RegularFixed<dimension>> mesh = std::dynamic_pointer_cast<Mesh::RegularFixed<dimension>>(data.mesh());
      
    // get data of solution vector
    std::vector<double> vectorValues;
    PetscUtility::getVectorEntries(data.solution(), vectorValues);
    
    // determine number of entries in each dimension
    std::vector<long int> nEntries(dimension);
    for (int i=0; i<dimension; i++)
    {
      nEntries[i] = (mesh->nElements(i) + 1) * data.nDegreesOfFreedomPerNode();
    }

    callCallback(vectorValues, nEntries, timeStepNo, currentTime);
  }
}

};