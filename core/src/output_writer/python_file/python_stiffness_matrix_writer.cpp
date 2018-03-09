#include "output_writer/python_file/python_stiffness_matrix_writer.h"

#include <Python.h>  // has to be the first included header
#include <memory>
#include <sstream>

namespace OutputWriter
{

void NumpyFileWriter::writeToNumpyFile(std::vector<double> &data, std::string filename, std::vector<long> &nEntries)
{
#if 1
  long int nEntriesTotal = 1;
  int dimension = nEntries.size();
  std::stringstream shape;
  shape << "[";
  for (int i=0; i<dimension; i++)
  {
    nEntriesTotal *= nEntries[i];
    if (i != 0)
      shape << ",";
    shape << nEntries[i];
    VLOG(1) << "nEntries["<<i<<"] = "<<nEntries[i];
  }
  shape << "]";
  
  if (nEntriesTotal != (long int)data.size())
  {
    LOG(ERROR) << "Number of entries " << nEntriesTotal << " " << nEntries << " for file \"" << filename << "\" does not match vector size "<<data.size()<<".";
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
  
  VLOG(1) << "data " << data << " written to file temp1";
  
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
    VLOG(1) << "converterScript: ["<<converterScript.str()<<"]";
    
    int ret = PyRun_SimpleString(converterScript.str().c_str());
    if (ret != 0)
      LOG(WARNING) << "Conversion to numpy file \"" << filename << "\" failed.";
    else
      LOG(INFO) << "Array of shape " << shape.str() << " exported to \"" << filename << "\"";
  }
  // remove temporary file
  //std::remove("temp1");
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
  
  
  PyObject *solutionVector = PyArray_SimpleNewFromData(nEntries.size(), nEntries.data(), NPY_DOUBLE, data.data());
  
  // write numpy array object to file
  PyObject *filenamePython = PyString_FromString(filename.c_str());
  PyArray_Dump(solutionVector, filenamePython, -1);
  
  Py_CLEAR(solutionVector);
  Py_CLEAR(filenamePython);
#endif

}

};   // namespace
