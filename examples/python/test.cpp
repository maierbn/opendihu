#include "opendihu.h"

#include <iostream>

int main(int argc, char* argv[])
{
  
  //wchar_t *programNameWChar = Py_DecodeLocale(argv[0], NULL);
  //Py_SetProgramName(programNameWChar);  /* optional but recommended */
    
  //std::string pythonSearchPath = std::string("/store/software/opendihu/dependencies/python/install");
  //const wchar_t *pythonSearchPathWChar = Py_DecodeLocale(pythonSearchPath.c_str(), NULL);
  //Py_SetPythonHome((wchar_t *)pythonSearchPathWChar);
  
  //std::string pythonPath("/store/software/opendihu/dependencies/python/install/lib/python3.6");
  //const wchar_t *pythonPathWChar = Py_DecodeLocale(pythonPath.c_str(), NULL);
  //Py_SetPath((wchar_t *)pythonPathWChar);
  
  std::cout<<"1) DihuContext settings"<<std::endl;
  DihuContext settings(argc, argv);
  
  std::cout<<"2) create data"<<std::endl;
  //Data::TimeStepping<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<1>,BasisFunction::LagrangeOfOrder<1>>, 1> data(settings.getPythonConfig());
  
  std::cout<<"3) test"<<std::endl;
  
  // test output
  PyObject *ioModule = PyImport_ImportModule("io");

  // open file
  std::string filename("out.py");
  PyObject *file = PyObject_CallMethod(ioModule, "open", "ss", filename.c_str(), "w");
  Py_DECREF(ioModule);
  
  // create test dict
  PyObject *pyData = PyDict_New();
  int ret = 0;
  ret = PyDict_SetItemString(pyData,"a", PyLong_FromLong(5));
  ret = PyDict_SetItemString(pyData,"b", PyUnicode_FromString("hi"));
  
  PythonUtility::printDict(pyData);
  
  // load json module
  static PyObject *jsonModule = NULL;
  if (jsonModule == NULL)
  {
    jsonModule = PyImport_ImportModule("json");
  }
  if (jsonModule == NULL)
  {
    LOG(ERROR) << "Could not import json module";
  }
  else
  {
    // convert data object to string representation and save in file
    //PyObject *pyDataJson = 
    PyObject_CallMethod(jsonModule, "dump", "O, O", pyData, file);
   
    if (PyErr_Occurred())
    {
      // print error message
      PyErr_Print();
    }
    
    PyObject_CallMethod(file, "write", "s", "W-ritten from Python C API!\n");
    
    
    if (PyErr_Occurred())
    {
      // print error message
      PyErr_Print();
    }
    
    PyObject_CallMethod(file, "flush", NULL);
    PyObject_CallMethod(file, "close", NULL);
    
    if (PyErr_Occurred())
    {
      // print error message
      PyErr_Print();
    }
  }
  
  Py_DECREF(file);
  // test output end
  //Py_Initialize();
  
  
  //Py_SetStandardStreamEncoding(NULL, NULL);
  
  //wchar_t **argvWChar = new wchar_t *[argc];
  //for (int i=0; i<argc; i++)
  //{
  //  argvWChar[i] = Py_DecodeLocale(argv[i], NULL);
  //}

  //PySys_SetArgvEx(argc, argvWChar, 1);
  //PyEval_InitThreads();

  std::cout<<"4) python version"<<std::endl;
  const char *version = Py_GetVersion();
  std::cout << "python version: " << version << std::endl;
    
  std::cout << std::string(80, '-') << std::endl;
  
  
  std::cout<<"5) load config"<<std::endl;
  std::string text = "print(\"hello\");config={\"format\": \"PythonFile\", \"filename\": \"out_diffusion2d\", \"binary\": False}";
  //int ret = 0;
  try
  {
    ret = PyRun_SimpleString(text.c_str());
  }
  catch(...)
  {
  }
  std::cout << std::string(80, '-') << std::endl;
  
  // if there was an error in the python code
  if (ret != 0)
  {
    if (PyErr_Occurred())
    {
      // print error message and exit
      PyErr_Print();
      exit(EXIT_FAILURE);
    }
    exit(EXIT_FAILURE);
  }
  std::cout<<"6) extract config"<<std::endl;
  
  PyObject *mainModule = PyImport_AddModule("__main__");
  PyObject *pythonConfig = PyObject_GetAttrString(mainModule, "config");
  
  std::cout<<"7) create pythonFile outputWriter"<<std::endl;
  OutputWriter::PythonFile pythonFile(pythonConfig); 
  
  
  std::cout<<"8) run write"<<std::endl;
  //pythonFile.write(data);
  std::cout<<"j"<<std::endl;
  
  
  
  
  
  
  
  
  
  
  std::string pythonConfigStr = R"(
# Diffusion 1D
n = 5
config = {
  "ExplicitEuler" : {
    "initialValues": [2,2,4,5,2,2],
    "numberTimeSteps": 5,
    "endTime": 1.0,
    "FiniteElementMethod" : {
      "nElements": n,
      "physicalExtent": 4.0,
      "relativeTolerance": 1e-15,
      "diffusionTensor": [5.0],
    },
    "OutputWriter" : [
      #{"format": "Paraview", "outputInterval": 1, "filename": "out", "binary": "false", "fixedFormat": False},
      {"format": "PythonFile", "filename": "out_diffusion1d", "outputInterval": 1, "binary":False}
    ]
  },
}
)";

  DihuContext settings2(argc, argv, pythonConfigStr);
  
  TimeSteppingScheme::ExplicitEuler<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredRegularFixedOfDimension<1>,
      BasisFunction::LagrangeOfOrder<>,
      Quadrature::None,
      Equation::Dynamic::IsotropicDiffusion
    >
  > problem(settings2);
  
  problem.run();
  
  
  
  
  

  return EXIT_SUCCESS;
}
