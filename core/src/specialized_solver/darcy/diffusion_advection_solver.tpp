#include "specialized_solver/darcy/diffusion_advection_solver.h"

#include <omp.h>
#include <sstream>

template<typename FiniteElementMethod>
DiffusionAdvectionSolver<FiniteElementMethod>::
DiffusionAdvectionSolver(DihuContext context) :
  Runnable(),
  ::TimeSteppingScheme::TimeSteppingScheme(context["DiffusionAdvectionSolver"]),   // replace "DiffusionAdvectionSolver" by the name of your solver, this will be the key for the dict in settings
  finiteElementMethod_(this->context_),
  data_(this->context_)
{
  // get python settings object from context
  this->specificSettings_ = this->context_.getPythonConfig();

  // initialize output writers
  this->outputWriterManager_.initialize(this->context_, this->specificSettings_);

  // parse options
  advectionVelocity_ = this->specificSettings_.getOptionArray<double,2>("advectionVelocity", std::array<double,2>{0,0});
  LOG(INFO) << "advectionVelocity: " << advectionVelocity_;
}

template<typename FiniteElementMethod>
void DiffusionAdvectionSolver<FiniteElementMethod>::
advanceTimeSpan(bool withOutputWritersEnabled)
{
  // This method computes some time steps of the simulation by running a for loop over the time steps.
  // The number of steps, timestep width and current time are all set by the parent class, TimeSteppingScheme.
  // You shouldn't change too much in this method.

  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  // compute time span of this method
  double timeSpan = this->endTime_ - this->startTime_;

  // output for debugging
  LOG(DEBUG) << "DiffusionAdvectionSolver::advanceTimeSpan, timeSpan=" << timeSpan<< ", timeStepWidth=" << this->timeStepWidth_
    << " n steps: " << this->numberTimeSteps_;

  Vec &solution = this->data_.solution()->getValuesContiguous();
  Vec &increment = this->data_.increment()->getValuesContiguous();

  // loop over time steps
  double currentTime = this->startTime_;
  for (int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    // in defined intervals (settings "timeStepOutputInterval") print out the current timestep
    if (timeStepNo % this->timeStepOutputInterval_ == 0 && (this->timeStepOutputInterval_ <= 10 || timeStepNo > 0))  // show first timestep only if timeStepOutputInterval is <= 10
    {
      LOG(INFO) << "DiffusionAdvectionSolver, timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime;
    }

    // Now, define what the solver does in this time step. Typically you want to execute the nested timestepping object.

    // Set timespan for finiteElementMethod_, the nested solver should advance its simulation by our timeStepWidth_.
    // This, in turn, may lead to multiple timesteps in the finiteElementMethod_.
    //this->finiteElementMethod_.setTimeSpan(currentTime, currentTime+this->timeStepWidth_);

    // advance the simulation by the specified time span, the parameter withOutputWritersEnabled specifies if the output writers should be called (usually yes)
    //finiteElementMethod_.advanceTimeSpan(withOutputWritersEnabled);

    // probably do something more here, maybe in a separate method:
    //executeMyHelperMethod();

    // ...

    Mat stiffnessMatrix = this->finiteElementMethod_.data().stiffnessMatrix()->valuesGlobal();
    //Mat massMatrix = this->finiteElementMethod_.data().massMatrix()->valuesGlobal();
    Mat inverseLumpedMassMatrix = this->finiteElementMethod_.data().inverseLumpedMassMatrix()->valuesGlobal();

    // convection term matrix
    Mat vMatrix = this->data_.vMatrix()->valuesGlobal();

    Mat KminusV;  // = K - V
    PetscErrorCode ierr;
    ierr = MatDuplicate(stiffnessMatrix, MAT_COPY_VALUES, &KminusV); CHKERRV(ierr);
    ierr = MatAXPY(KminusV, -1, vMatrix, SAME_NONZERO_PATTERN); CHKERRV(ierr);

    // compute M^-1*(K - V)
    Mat inverseMKminusV;

    ierr = MatMatMult(inverseLumpedMassMatrix, KminusV, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &inverseMKminusV); CHKERRV(ierr);

    // increment = (M^-1*(K - V)) * u
    ierr = MatMult(inverseMKminusV, solution, increment); CHKERRV(ierr);

    // u = u + dt * increment
    // y = alpha x + y.
    ierr = VecAXPY(solution, this->timeStepWidth_, increment); CHKERRV(ierr);


    // advance simulation time
    timeStepNo++;

    // compute new current simulation time
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;

    // stop duration measurement
    if (this->durationLogKey_ != "")
      Control::PerformanceMeasurement::stop(this->durationLogKey_);

    // if the output writers are enabled, write current output values using the output writers
    if (withOutputWritersEnabled)
      this->outputWriterManager_.writeOutput(this->data_, timeStepNo, currentTime);

    // start duration measurement
    if (this->durationLogKey_ != "")
      Control::PerformanceMeasurement::start(this->durationLogKey_);
  }

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);
}

template<typename FiniteElementMethod>
void DiffusionAdvectionSolver<FiniteElementMethod>::
initialize()
{
  // initialize() will be called before the simulation starts.

  // call initialize of the parent class, this parses the timestepping settings from the settings file
  TimeSteppingScheme::TimeSteppingScheme::initialize();

  // add this solver to the solvers diagram, which is an ASCII art representation that will be created at the end of the simulation.
  DihuContext::solverStructureVisualizer()->addSolver("DiffusionAdvectionSolver", true);   // hasInternalConnectionToFirstNestedSolver=true (the last argument) means slot connector data is shared with the first subsolver
  // if you have your own slot connector data rather than the one of the subsolver, call "addSolver" with false as second argument

  // indicate in solverStructureVisualizer that now a child solver will be initialized
  DihuContext::solverStructureVisualizer()->beginChild();

  // call initialize of the nested timestepping solver
  finiteElementMethod_.initialize();
  finiteElementMethod_.initializeForImplicitTimeStepping();

  // indicate in solverStructureVisualizer that the child solver initialization is done
  DihuContext::solverStructureVisualizer()->endChild();

  // In order to initialize the data object and actuall create all variables, we first need to assign a function space to the data object.
  // A function space object of type FunctionSpace<MeshType,BasisFunctionType> (see "function_space/function_space.h")
  // is an object that stores the mesh (e.g., all nodes and elements) as well as the basis function (e.g. linear Lagrange basis functions).
  // The finiteElementMethod_ solver already created a function space that we should use. We already have a typedef "FunctionSpace" that is the class of finiteElementMethod_'s function space type.
  std::shared_ptr<FunctionSpace> functionSpace = finiteElementMethod_.data().functionSpace();

  // Pass the function space to the data object. data_ stores field variables.
  // It needs to know the number of nodes and degrees of freedom (dof) from the function space in order to create the vectors with the right size.
  data_.setFunctionSpace(functionSpace);

  // now call initialize, data will then create all variables (Petsc Vec's)
  data_.initialize();

  // set the slotConnectorData for the solverStructureVisualizer to appear in the solver diagram
  DihuContext::solverStructureVisualizer()->setSlotConnectorData(getSlotConnectorData());

  // initialize values of solution vector

  // get number of global dofs, i.e. number of values in global list
  const int nDofsGlobal = this->data_.functionSpace()->nDofsGlobal();

  // set initial values as given in settings, or set to zero if not given
  std::vector<VecD<1>> localValues;
  this->specificSettings_.getOptionVector("initialValues", nDofsGlobal, localValues);

  // extract only the local dofs out of the list of global values
  this->data_.functionSpace()->meshPartition()->extractLocalDofsWithoutGhosts(localValues);

  this->data_.solution()->setValuesWithoutGhosts(localValues);

  // compute V system matrix
  computeVMatrix();
}

template<typename FiniteElementMethod>
void DiffusionAdvectionSolver<FiniteElementMethod>::
run()
{
  // The run method should not be changed. It is the method that gets called directly from the example main file.
  // If this solver itself is nested in other solvers or coupling schemes,
  // run() will not be called, but the enclosing solver will call initialize() and advanceTimeSpan().
  initialize();

  advanceTimeSpan();
}


template<typename FiniteElementMethod>
void DiffusionAdvectionSolver<FiniteElementMethod>::
computeVMatrix()
{

  typedef Quadrature::Gauss<2> QuadratureType;
  const int nComponents = 1;

  const int D = FunctionSpace::dim();
  LOG(TRACE) << "setStiffnessMatrix " << D << "D using integration, FunctionSpace: " << StringUtility::demangle(typeid(FunctionSpace).name())
    << ", QuadratureType: " << StringUtility::demangle(typeid(QuadratureType).name());

  // define shortcuts for integrator and basis
  typedef Quadrature::TensorProduct<D,QuadratureType> QuadratureDD;
  const int nDofsPerElement = FunctionSpace::nDofsPerElement();
  const int nUnknownsPerElement = nDofsPerElement*nComponents;
  typedef MathUtility::Matrix<nUnknownsPerElement,nUnknownsPerElement,double_v_t> EvaluationsType;
  typedef std::array<
            EvaluationsType,
            QuadratureDD::numberEvaluations()
          > EvaluationsArrayType;     // evaluations[nGP^D][nDofs][nDofs]

  // setup arrays used for integration
  std::array<std::array<double,D>, QuadratureDD::numberEvaluations()> samplingPoints = QuadratureDD::samplingPoints();
  EvaluationsArrayType evaluationsArray{};

  LOG(DEBUG) << "1D integration with " << QuadratureType::numberEvaluations() << " evaluations";
  LOG(DEBUG) << D << "D integration with " << QuadratureDD::numberEvaluations() << " evaluations";

  // initialize variables
  std::shared_ptr<PartitionedPetscMat<FunctionSpace>> vMatrix = this->data_.vMatrix();

  std::shared_ptr<FunctionSpace> functionSpace = std::static_pointer_cast<FunctionSpace>(this->data_.functionSpace());
  functionSpace->geometryField().setRepresentationGlobal();
  functionSpace->geometryField().startGhostManipulation();   // ensure that local ghost values of geometry field are set

  bool outputAssemble3DStiffnessMatrixHere = false;
  std::chrono::time_point<std::chrono::system_clock> tStart;

  const element_no_t nElementsLocal = functionSpace->nElementsLocal();
  LOG(DEBUG) << " nElementsLocal: " << nElementsLocal;

  // initialize values to zero
  // loop over elements, always 4 elements at once using the vectorized functions
  for (int elementNoLocal = 0; elementNoLocal < nElementsLocal; elementNoLocal += nVcComponents)
  {

#ifdef USE_VECTORIZED_FE_MATRIX_ASSEMBLY
    // get indices of elementNos that should be handled in the current iterations,
    // this is, e.g.
    //    [10,11,12,13,-1,-1,-1,-1] (if nVcComponents==4 and nElementsLocal > 13)
    // or [10,11,12,-1,-1,-1,-1,-1] (if nVcComponents==4 and nElementsLocal == 13)

    dof_no_v_t elementNoLocalv([elementNoLocal, nElementsLocal](int i)
    {
      return (i >= nVcComponents || elementNoLocal+i >= nElementsLocal? -1: elementNoLocal+i);
    });

    // here, elementNoLocalv is the list of indices of the current iteration, e.g. [10,11,12,13,-1,-1,-1,-1]
    // elementNoLocal is the first entry of elementNoLocalv
#else
    int elementNoLocalv = elementNoLocal;
#endif

    std::array<dof_no_v_t,nDofsPerElement> dofNosLocal = functionSpace->getElementDofNosLocal(elementNoLocalv);

    for (int i = 0; i < nDofsPerElement; i++)
    {
      for (int j = 0; j < nDofsPerElement; j++)
      {
        // loop over components (1,...,D for solid mechanics)
        for (int rowComponentNo = 0; rowComponentNo < nComponents; rowComponentNo++)
        {
          for (int columnComponentNo = 0; columnComponentNo < nComponents; columnComponentNo++)
          {
            int componentNo = rowComponentNo*nComponents + columnComponentNo;

            //LOG(DEBUG) << " initialize vMatrix entry ( " << dofNosLocal[i] << "," << dofNosLocal[j] << ") (no. " << cntr++ << ")";
            vMatrix->setValue(componentNo, dofNosLocal[i], dofNosLocal[j], 0, INSERT_VALUES);
          }
        }
      }
    }
  }

  // allow switching between vMatrix->setValue(... INSERT_VALUES) and ADD_VALUES
  vMatrix->assembly(MAT_FLUSH_ASSEMBLY);

  double progress = 0;

  // fill entries in stiffness matrix
  // loop over elements, always 4 elements at once using the vectorized functions
  for (int elementNoLocal = 0; elementNoLocal < nElementsLocal; elementNoLocal += nVcComponents)
  {

#ifdef USE_VECTORIZED_FE_MATRIX_ASSEMBLY
    // get indices of elementNos that should be handled in the current iterations,
    // this is, e.g.
    //    [10,11,12,13,-1,-1,-1,-1] (if nVcComponents==4 and nElementsLocal > 13)
    // or [10,11,12,-1,-1,-1,-1,-1] (if nVcComponents==4 and nElementsLocal == 13)

    dof_no_v_t elementNoLocalv([elementNoLocal, nElementsLocal](int i)
    {
      return (i >= nVcComponents || elementNoLocal+i >= nElementsLocal? -1: elementNoLocal+i);
    });

    // here, elementNoLocalv is the list of indices of the current iteration, e.g. [10,11,12,13,-1,-1,-1,-1]
    // elementNoLocal is the first entry of elementNoLocalv
#else
    int elementNoLocalv = elementNoLocal;
#endif

    if (outputAssemble3DStiffnessMatrixHere && this->context_.ownRankNoCommWorld() == 0)
    {
      double newProgress = (double)elementNoLocal / nElementsLocal;
      if (int(newProgress*10) != int(progress*10))
      {
        std::cout << "\b\b\b\b" << int(newProgress*100) << "%" << std::flush;
      }
      progress = newProgress;
    }

    // get indices of element-local dofs
    std::array<dof_no_v_t,nDofsPerElement> dofNosLocal = functionSpace->getElementDofNosLocal(elementNoLocalv);

    VLOG(2) << "element " << elementNoLocalv;

    // get geometry field (which are the node positions for Lagrange basis and node positions and derivatives for Hermite)
    std::array<Vec3_v_t,FunctionSpace::nDofsPerElement()> geometry;
    functionSpace->getElementGeometry(elementNoLocalv, geometry);

    // compute integral
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPoints.size(); samplingPointIndex++)
    {
      // evaluate function to integrate at samplingPoint
      std::array<double,D> xi = samplingPoints[samplingPointIndex];

      // compute the 3xD jacobian of the parameter space to world space mapping
      std::array<Vec3_v_t,D> jacobian = FunctionSpace::computeJacobian(geometry, xi);

      // get the factor in the integral that arises from the change in integration domain from world to coordinate space
      double_v_t integrationFactor = MathUtility::computeIntegrationFactor(jacobian);

      VLOG(2) << "samplingPointIndex=" << samplingPointIndex<< ", xi=" <<xi<< ", geometry: " <<geometry<< ", jac: " <<jacobian;

      std::array<std::array<double,D>,FunctionSpace::nDofsPerElement()> gradPhi = data_.functionSpace()->getGradPhi(xi);

      // get evaluations of integrand at xi for all (i,j)-dof pairs, integrand is defined in another class
      for (int i = 0; i < nDofsPerElement; i++)
      {
        for (int j = 0; j < nDofsPerElement; j++)
        {
          // v * phi[i] * gradPhi[j](xi)
          evaluationsArray[samplingPointIndex](i,j)
            = FunctionSpace::phi(i,xi) * (gradPhi[j][0] * advectionVelocity_[0] + gradPhi[j][1] * advectionVelocity_[1]) * integrationFactor;
        }
      }
    }  // function evaluations

    // integrate all values for the (i,j) dof pairs at once
    EvaluationsType integratedValues = QuadratureDD::computeIntegral(evaluationsArray);

    // perform integration and add to entry of stiffness matrix
    for (int i = 0; i < nDofsPerElement; i++)
    {
      for (int j = 0; j < nDofsPerElement; j++)
      {
        // loop over components (1,...,D for solid mechanics)
        for (int rowComponentNo = 0; rowComponentNo < nComponents; rowComponentNo++)
        {
          for (int columnComponentNo = 0; columnComponentNo < nComponents; columnComponentNo++)
          {
            // integrate value and set entry in stiffness matrix
            double_v_t integratedValue = integratedValues(i*nComponents + rowComponentNo, j*nComponents + columnComponentNo);
            int componentNo = rowComponentNo*nComponents + columnComponentNo;

            VLOG(2) << "  dof pair (" << i<< "," <<j<< ") dofs (" << dofNosLocal[i]<< "," << dofNosLocal[j]<< "), "
              << "component (" << rowComponentNo << "," << columnComponentNo << "), " << componentNo
              << ", integrated value: " <<integratedValue;

            // get local dof no
            dof_no_v_t dofINoLocal = dofNosLocal[i];
            dof_no_v_t dofJNoLocal = dofNosLocal[j];

            // add the entry in the stiffness matrix, for all dofs of the vectorized values at once,
            // i.e. K_dofINoLocal[0],dofJNoLocal[0] = value[0]
            // i.e. K_dofINoLocal[1],dofJNoLocal[1] = value[1], etc.
            // Note that K_dofINoLocal[0],dofJNoLocal[1] would be potentially zero, the contributions are considered element-wise
            vMatrix->setValue(componentNo, dofINoLocal, dofJNoLocal, integratedValue, ADD_VALUES);
          }
        }
      }  // j
    }  // i
  }  // elementNoLocalv

  vMatrix->assembly(MAT_FINAL_ASSEMBLY);
}

template<typename FiniteElementMethod>
void DiffusionAdvectionSolver<FiniteElementMethod>::
reset()
{
  finiteElementMethod_.reset();

  // "uninitialize" everything
}

//! call the output writer on the data object, output files will contain currentTime, with callCountIncrement !=1 output timesteps can be skipped
template<typename FiniteElementMethod>
void DiffusionAdvectionSolver<FiniteElementMethod>::
callOutputWriter(int timeStepNo, double currentTime, int callCountIncrement)
{
  // Call the output writer manager to execute all output writers.
  // This will check if there is a file to be written at the current timestep and then output the files.
  this->outputWriterManager_.writeOutput(this->data_, timeStepNo, currentTime, callCountIncrement);
}

template<typename FiniteElementMethod>
void DiffusionAdvectionSolver<FiniteElementMethod>::
executeMyHelperMethod()
{
  // this is the template for an own private method

  // for example you can get the solution values of the finiteElementMethod_ by
  Vec solution = finiteElementMethod_.data().solution()->valuesGlobal();

  // As an example we invert all the solution values.
  // Because "Vec"'s are actually pointers, this effects the actual data, no copy-back is needed.
  // Note the error handling with Petsc functions. Always use "ierr = PetscFunction(); CHKERRV(ierr);"
  PetscErrorCode ierr;
  ierr = VecScale(solution, -1); CHKERRV(ierr);

  // get a field variable from data object:
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace,1>> fieldVariableA = data_.fieldVariableA();

  // copy the solution from the finiteElementMethod_ to fieldVariableA.
  // note, you get the Petsc "Vec" of a field variable by "valuesGlobal()"
  ierr = VecCopy(solution, fieldVariableA->valuesGlobal()); CHKERRV(ierr);

  // add 1.0 to every entry in fieldVariableA, this already updates fieldVariableA in data because it is a pointer. There is no need to copy the values back.
  ierr = VecShift(fieldVariableA->valuesGlobal(), 1.0); CHKERRV(ierr);
}

template<typename FiniteElementMethod>
typename DiffusionAdvectionSolver<FiniteElementMethod>::Data &DiffusionAdvectionSolver<FiniteElementMethod>::
data()
{
  // get a reference to the data object
  return data_;

  // The finiteElementMethod_ object also has a data object, we could also directly use this and avoid having an own data object:
  //  return finiteElementMethod_.data();
}

//! get the data that will be transferred in the operator splitting to the other term of the splitting
//! the transfer is done by the slot_connector_data_transfer class
template<typename FiniteElementMethod>
std::shared_ptr<typename DiffusionAdvectionSolver<FiniteElementMethod>::SlotConnectorDataType> DiffusionAdvectionSolver<FiniteElementMethod>::
getSlotConnectorData()
{
  //! This is relevant only, if this solver is part of a splitting or coupling scheme. Then this method returns the values/variables that will be
  // transferred to the other solvers. We can just reuse the values of the finiteElementMethod_.
  return finiteElementMethod_.getSlotConnectorData();
}
