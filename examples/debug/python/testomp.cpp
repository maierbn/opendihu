#include "opendihu.h"

#include <iostream>
#include <sstream>

#include <omp.h>

int main(int argc, char* argv[])
{
  std::string pythonConfig = R"(
config = {
  "Solvers": {
    "testname": {
      "relativeTolerance": 1e-15,
    }
  },
  "ExplicitEuler" : {
    "initialValues": [2,2,4,5,2,2],
    "numberTimeSteps": 5000,
    "endTime": 1.0,
    "FiniteElementMethod" : {
      "nElements": 4,
      "solverName": "testname",
      "physicalExtent": 4.0,
      "initialValues": [0],
      "DirichletBoundaryCondition": {0:1.0},
      "nodePositions": [[0,0,0], [1,0], [2,0,0], [0,1], [1,1], [2,1], [0,2], [1,2], [2,2]],
      "elements": [[0, 1, 3, 4], [1, 2, 4, 5], [3, 4, 6, 7], [4, 5, 7, 8]],  
    },
  }
}
)";
  DihuContext context(argc, argv, pythonConfig);
  /*
  TimeSteppingScheme::ExplicitEuler<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredRegularFixedOfDimension<1>,
      BasisFunction::LagrangeOfOrder<>,
      Quadrature::None,
      Equation::Dynamic::IsotropicDiffusion
    >
  > problem(context);
  problem.run();
  */
  /*
  Control::MultipleInstances<
    TimeSteppingScheme::ExplicitEuler<
      SpatialDiscretization::FiniteElementMethod<
        Mesh::StructuredRegularFixedOfDimension<1>,
        BasisFunction::LagrangeOfOrder<1>,
        Quadrature::Gauss<2>,
        Equation::Dynamic::IsotropicDiffusion
      >
    >
  > problem(context);
  
  problem.run();
  */
  
  typedef TimeSteppingScheme::ExplicitEuler<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredRegularFixedOfDimension<1>,
      BasisFunction::LagrangeOfOrder<1>,
      Quadrature::Gauss<2>,
      Equation::Dynamic::IsotropicDiffusion
    >
  > ProblemType;
  
  std::vector<ProblemType> instances_(2,context);
  
  //instances_[0].initialize();
  //instances_[1].initialize(); 
  //instances_[0].advanceTimeSpan();
  
  
  //problem.initialize();
 
  int nInstances_ = 2;
  #pragma omp parallel for
  for (int i = 0; i < nInstances_; i++)
  {
    if (omp_get_thread_num() == 0)
    {
      std::stringstream msg;
      msg << omp_get_thread_num() << ": running " << nInstances_ << " instances with " << omp_get_num_threads() << " OpenMP threads";
      LOG(INFO) << msg.str();
    }
    
    
    
    if (!instances_[i].initialized_)
    {
      
      PyGILState_STATE gstate_;
      gstate_ = PyGILState_Ensure();
      
      LOG(TRACE) << omp_get_thread_num() << ": TimeSteppingSchemeOde::initialize -*";

      // initialize underlying DiscretizableInTime object
      instances_[i].discretizableInTime_.initialize();

      LOG(TRACE) << omp_get_thread_num() << ": set mesh";
      std::shared_ptr<Mesh::Mesh> mesh = instances_[i].discretizableInTime_.mesh();
      instances_[i].data_->setMesh(std::static_pointer_cast<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<1>, BasisFunction::LagrangeOfOrder<1>>>(mesh));
      instances_[i].data_->initialize();

      instances_[i].timeStepOutputInterval_ = PythonUtility::getOptionInt(instances_[i].specificSettings_, "timeStepOutputInterval", 100, PythonUtility::Positive);

      // set initial values from settings

      Vec &solution = instances_[i].data_->solution().values();

      if (!instances_[i].discretizableInTime_.setInitialValues(solution))
      {
        instances_[i].setInitialValues();
      }
      
      instances_[i].initialized_ = true;
      
    
      PyGILState_Release(gstate_);
      
      
      
    }
      
    instances_[i].run();
  }
  
  
  exit(0);
  // test output
  const char *version = Py_GetVersion();
  std::cout << "python version: " << version << std::endl;
   
  #pragma omp parallel for
  for (int i=0; i<2; i++)
  {
    if (i == 0)
    {
      double k = 0;
      for (int j=0; j<10000000; j++)
      {
        k += sqrt(j);
      }
      std::cout<<" k = " <<k<<std::endl;
    }
    
    PyObject *settings = context.getPythonConfig();
    
    for (int j=0; j<10000; j++)
    {
      std::cout << int(omp_get_thread_num()) << "/" << int(omp_get_num_threads()) << " start i=" << i << ", j=" << j << std::endl;
      
      PyGILState_STATE gstate_;
      gstate_ = PyGILState_Ensure();
              
      
              
      LOG(INFO) << omp_get_thread_num() << ": discretizableInTime::initializeLinearSolver";
    
      //problem.discretizableInTime().initializeLinearSolver();
      
      // if solver has already been created earlier
      if (PythonUtility::hasKey(settings, "solverName"))
      {
        LOG(INFO) << omp_get_thread_num() << ": Manager::solver - hasKey";
      }
      
    
      PyObject *pyData = PyDict_New();
      PyDict_SetItemString(pyData, "a", PyLong_FromLong(5));
      PyDict_SetItemString(pyData, "b", PyUnicode_FromString("hi"));
      
      PyObject *ioModule = PyImport_ImportModule("io");

      PyGILState_STATE gstate2_;
      gstate2_ = PyGILState_Ensure();
      
      std::stringstream filename;
      filename << "out/testio_" << i << "_" << j << ".py";
      PyObject *file = PyObject_CallMethod(ioModule, "open", "ss", filename.str().c_str(), "w");
      Py_DECREF(ioModule);
      
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
        
        PyObject_CallMethod(file, "flush", NULL);
        PyObject_CallMethod(file, "close", NULL);
      }
      
      
      PyGILState_Release(gstate2_);
      
      PyGILState_Release(gstate_);
    }
  }
  
  std::cout << "done." << std::endl;
  
  return EXIT_SUCCESS;
}
