#include "output_writer/python.h"

#include <fstream>
#include <iostream>
#include <thread>
#include <chrono>

#include "easylogging++.h"
#include <Python.h>
#include <numpy/ndarraytypes.h>
#include <numpy/npy_common.h>
#include <numpy/npy_math.h>
#include <numpy/ndarrayobject.h>

#include <control/python_utility.h>
#include <mesh/regular_fixed.h>
#include <mesh/rectilinear_fixed.h>
#include <mesh/nonrectilinear_fixed.h>
#include <mesh/deformable.h>
#include <mesh/mesh.h>

namespace OutputWriter
{

Python::Python(PyObject *settings) : settings_(settings)
{

}

void Python::writeToNumpyFile(std::vector<double> &data, std::string filename, int dimension, std::vector<long> &nEntries)
{
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
  for(auto &value : data)
  {
    union {
      double d;
      char c[8];
    };
    d = value;
    for(int i=0; i<8; i++)
    file<<c[i];
  }
  
  file.close();
  
  // convert to numpy file by python script
  std::stringstream converterScript;
  converterScript << "import numpy as np" << std::endl
    << "v = np.fromfile(\"temp1\")" << std::endl 
    << "v = np.reshape(v," << shape.str() << ")" << std::endl
    << "np.save(\"" << filename << "\",v)";
  
  int ret = PyRun_SimpleString(converterScript.str().c_str());
  if (ret != 0)
    LOG(WARNING) << "Conversion to numpy file \"" << filename << "\" failed.";
  else
    LOG(INFO) << "Array of shape " << shape.str() << " exported to \"" << filename << "\"";
  
  // remove temporary file
  std::remove("temp1");
    
  // directly write npy file by using numpy c API (not working)
#if 0    
  // construct numpy array object
  long int nEntriesTotal = 1;
  std::vector<long int> nEntries(dimension);
  for (int i=0; i<dimension; i++)
  {
    nEntries[i] = (mesh->nElements(i) + 1);
    nEntriesTotal *= nEntries[i];
    LOG(DEBUG) << "nEntries["<<i<<"] = "<<nEntries[i];
  }
  
  if (nEntriesTotal != vectorSize) 
  {
    LOG(ERROR) << "number of degrees of freedom " << nEntriesTotal << " does not match vector size of solution "<<vectorSize<<".";
    return;
  }
  
  LOG(DEBUG) << "vectorSize="<<vectorSize<<", "<<vectorValues.size()<<", nEntries size "<<nEntries.size();
  
  // test PyArray_SimpleNewFromData
  long int dims[2] = {2,2};
  long int ds[1] = {1};
  double array[2][2] = {1.0, 2.0, 3.0, 4.0};
  //PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, array);
  PyArray_SimpleNew(1, ds, NPY_INT);

  LOG(DEBUG) << "success";
  //PyObject *filename0 = PyString_FromString("test.npy");
  //PyArray_Dump(test, filename0, -1);
  
  
  PyObject *solutionVector = PyArray_SimpleNewFromData(dimension, nEntries.data(), NPY_DOUBLE, vectorValues.data());
  
  // write numpy array object to file
  PyObject *filename = PyString_FromString(filenameSolution.c_str());
  PyArray_Dump(solutionVector, filename, -1);
  
  Py_XDECREF(solutionVector);
  Py_XDECREF(filename);
#endif

}

void Python::writeSolution(Data::Data& data)
{
  const int dimension = data.mesh()->dimension();
  
  LOG(DEBUG) << "dimension: "<<dimension;
  
  // solution and rhs vectors in mesh shape
  switch(dimension)
  {
  case 1:
    writeSolutionDim<1>(data);
    break;
  case 2:
    writeSolutionDim<2>(data);
    break;
  case 3:
    writeSolutionDim<3>(data);
    break;
  };
}
template <int dimension>
void Python::writeSolutionDim(Data::Data &data)
{
  LOG(DEBUG) << "writeSolution<"<<dimension<<">()";
  
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
    int vectorSize = 0;
    VecGetSize(data.solution(), &vectorSize);

    std::vector<int> indices(vectorSize);
    std::iota(indices.begin(), indices.end(), 0);
    std::vector<double> vectorValues(vectorSize);

    VecGetValues(data.solution(), vectorSize, indices.data(), vectorValues.data());
    
    // determine number of entries in each dimension
    LOG(DEBUG) << "dimension=" << dimension;
    
    std::vector<long int> nEntries(dimension);
    for (int i=0; i<dimension; i++)
    {
      nEntries[i] = (mesh->nElements(i) + 1);
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
  LOG(DEBUG) << "writeRhsMatrix<"<<dimension<<">()";
  
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
    LOG(DEBUG) << "dimension=" << dimension;
    
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