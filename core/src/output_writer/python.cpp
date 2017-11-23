#include "output_writer/python.h"

#include <fstream>
#include <iostream>
#include <thread>
#include <chrono>

#include "easylogging++.h"
#include <Python.h>

// numpy api
#if 0
#include <numpy/ndarraytypes.h>
#include <numpy/npy_common.h>
#include <numpy/npy_math.h>
#include <numpy/ndarrayobject.h>
#endif

#include <control/python_utility.h>
#include <control/petsc_utility.h>
#include <mesh/regular_fixed.h>
#include <mesh/rectilinear_fixed.h>
#include <mesh/nonrectilinear_fixed.h>
#include <mesh/deformable.h>
#include <mesh/mesh.h>

namespace OutputWriter
{

Python::Python(PyObject *settings) : Generic(settings)
{

}

void Python::writeSolution(Data::Data& data, int timeStepNo, double currentTime)
{
  if (!data.mesh())
  {
    LOG(FATAL) << "mesh is not set!";
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

void Python::writeToNumpyFile(std::vector<double> &data, std::string filename, int dimension, std::vector<long> &nEntries)
{
#if 1
  long int nEntriesTotal = 1;
  std::stringstream shape;
  shape << "[";
  for (int i=0; i<dimension; i++)
  {
    nEntriesTotal *= nEntries[i];
    if (i != 0)
      shape << ",";
    shape << nEntries[i];
    LOG(DEBUG) << "nEntries["<<i<<"] = "<<nEntries[i];
  }
  shape << "]";
  
  if (nEntriesTotal != (long int)data.size())
  {
    LOG(ERROR) << "Number of entries " << nEntriesTotal << " does not match vector size "<<data.size()<<".";
    return;
  }
  
  // write data to binary file
  std::ofstream file("temp1", std::ios::out | std::ios::binary | std::ios::trunc);
  if (!file.is_open())
  {
    LOG(ERROR) << "Could not write temporary files!";
    return;
  }
  for(auto value : data)
  {
    union {
      double d;
      char c[8];
    };
    
    if (std::isnan(d))
      d = 0.0;
    else
      d = value;
    
    
    for(int i=0; i<8; i++)
    {
      file<<c[i];
    }
  }
  
  file.close();
  
  // convert to numpy file by python script
  std::stringstream converterScript;
  converterScript << "import numpy as np" << std::endl
    << "v = np.fromfile(\"temp1\")" << std::endl 
    << "v = np.reshape(v," << shape.str() << ")" << std::endl
    << "np.save(\"" << filename << "\",v)";
  
  //LOG(DEBUG) << converterScript.str();
  if (0)
  {
    std::ofstream scriptFile("convert.py");
    scriptFile<<converterScript.str();
    scriptFile.close();
    int ret = system("python convert.py");
    if (ret)
      LOG(DEBUG) << "convert script failed!";
    std::remove("convert.py");
  }
  else
  {
    int ret = PyRun_SimpleString(converterScript.str().c_str());
    if (ret != 0)
      LOG(WARNING) << "Conversion to numpy file \"" << filename << "\" failed.";
    else
      LOG(INFO) << "Array of shape " << shape.str() << " exported to \"" << filename << "\"";
  }
  // remove temporary file
  std::remove("temp1");
#endif    
  // directly write npy file by using numpy c API (not working)
#if 0
  // construct numpy array object
  long int nEntriesTotal = 1;
  for (int i=0; i<dimension; i++)
  {
    nEntriesTotal *= nEntries[i];
    LOG(DEBUG) << "nEntries["<<i<<"] = "<<nEntries[i];
  }
  
  if (nEntriesTotal != data.size()) 
  {
    LOG(ERROR) << "Number of entries " << nEntriesTotal << " does not match vector size "<<data.size()<<".";
    return;
  }
  
  // test PyArray_SimpleNewFromData
  long int dims[2] = {2,2};
  long int ds[1] = {1};
  double array[2][2] = {1.0, 2.0, 3.0, 4.0};
  //PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, array);
  PyArray_SimpleNew(1, ds, NPY_INT);

  LOG(DEBUG) << "success";
  //PyObject *filename0 = PyString_FromString("test.npy");
  //PyArray_Dump(test, filename0, -1);
  
  
  PyObject *solutionVector = PyArray_SimpleNewFromData(dimension, nEntries.data(), NPY_DOUBLE, data.data());
  
  // write numpy array object to file
  PyObject *filenamePython = PyString_FromString(filename.c_str());
  PyArray_Dump(solutionVector, filenamePython, -1);
  
  Py_XDECREF(solutionVector);
  Py_XDECREF(filenamePython);
#endif

}

template <int dimension>
void Python::writeSolutionDim(Data::Data &data, int timeStepNo, double currentTime)
{
  LOG(TRACE) << "writeSolution<"<<dimension<<">()";
  
  // solution and rhs vectors in mesh shape
  // if mesh is regular fixed
  if (std::dynamic_pointer_cast<Mesh::RegularFixed<dimension>>(data.mesh()) != NULL)
  {
    // cast mesh to its special type
    std::shared_ptr<Mesh::RegularFixed<dimension>> mesh = std::dynamic_pointer_cast<Mesh::RegularFixed<dimension>>(data.mesh());
      
    // determine file names
    std::stringstream s[2];
    s[0] << filename_ << "_solution.npy";
    s[1] << filename_ << "_solution_shaped.npy";
    std::string filenameSolution = s[0].str();
    std::string filenameSolutionShaped = s[1].str();
    
    // get data of solution vector
    std::vector<double> vectorValues;
    PetscUtility::getVectorEntries(data.solution(), vectorValues);
    
    // determine number of entries in each dimension
    std::vector<long int> nEntries(dimension);
    for (int i=0; i<dimension; i++)
    {
      nEntries[i] = (mesh->nElements(i) + 1) * data.nDegreesOfFreedomPerNode();
    }
    std::vector<long int> singleEntry({(long)vectorValues.size()});
    
    // write as numpy file
    writeToNumpyFile(vectorValues, filenameSolution, 1, singleEntry);
    writeToNumpyFile(vectorValues, filenameSolutionShaped, dimension, nEntries);
    
    // if data is of class Data::FiniteElements
    if (dynamic_cast<Data::FiniteElements *>(&data) != NULL)
    {
      writeRhsMatrix<dimension>(*dynamic_cast<Data::FiniteElements *>(&data));
    }
  }
}

template <int dimension>
void Python::writeRhsMatrix(Data::FiniteElements &data)
{
  // solution and rhs vectors in mesh shape
  if (std::dynamic_pointer_cast<Mesh::RegularFixed<dimension>>(data.mesh()) != NULL)
  {
    // cast mesh to its special type
    std::shared_ptr<Mesh::RegularFixed<dimension>> mesh = std::dynamic_pointer_cast<Mesh::RegularFixed<dimension>>(data.mesh());
      
    // determine file names
    std::stringstream s[3];
    s[0] << filename_ << "_rhs.npy";
    s[1] << filename_ << "_rhs_shaped.npy";
    s[2] << filename_ << "_stiffness.npy";
    std::string filenameRhs = s[0].str();
    std::string filenameRhsShaped = s[1].str();
    std::string filenameStiffness = s[2].str();
    
    // get data of rhs vector
    int vectorSize = 0;
    VecGetSize(data.rightHandSide(), &vectorSize);

    std::vector<int> indices(vectorSize);
    std::iota(indices.begin(), indices.end(), 0);
    std::vector<double> vectorValues(vectorSize);

    // determine number of entries in each dimension
    std::vector<long int> nEntries(dimension);
    for (int i=0; i<dimension; i++)
    {
      nEntries[i] = (mesh->nElements(i) + 1);
    }
    std::vector<long int> singleEntry({(long)vectorValues.size()});
    
    VecGetValues(data.rightHandSide(), vectorSize, indices.data(), vectorValues.data());
    
    // write as numpy file
    writeToNumpyFile(vectorValues, filenameRhs, 1, singleEntry);
    writeToNumpyFile(vectorValues, filenameRhsShaped, dimension, nEntries);
    
    // get stiffness matrix
    int nRows, nColumns;
    MatGetSize(data.stiffnessMatrix(), &nRows, &nColumns);
    std::vector<int> rowIndices(nRows);
    std::iota(rowIndices.begin(), rowIndices.end(), 0);
    std::vector<int> columnIndices(nColumns);
    std::iota(columnIndices.begin(), columnIndices.end(), 0);
    std::vector<double> matrixValues(nRows*nColumns);
    
    nEntries = {nRows, nColumns};
    
    MatGetValues(data.stiffnessMatrix(), nRows, rowIndices.data(), nColumns, columnIndices.data(), matrixValues.data());
    
    // write as numpy file
    writeToNumpyFile(matrixValues, filenameStiffness, 2, nEntries);
    
  }
}
};