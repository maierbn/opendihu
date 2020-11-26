#include "specialized_solver/solid_mechanics/quasi_static/quasi_static_nonlinear_elasticity_solver_chaste.h"

#include <Python.h>  // has to be the first included header

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "data_management/specialized_solver/multidomain.h"
#include "control/diagnostic_tool/performance_measurement.h"

#ifdef HAVE_CHASTE

#include "ExecutableSupport.hpp"
#include "Exception.hpp"
#include "PetscTools.hpp"
#include "PetscException.hpp"

#include <vector>

#include "UblasCustomFunctions.hpp"
/* The incompressible solver is called `IncompressibleNonlinearElasticitySolver` */
#include "IncompressibleNonlinearElasticitySolver.hpp"
/* The simplest incompressible material law is the Mooney-Rivlin material law (of which
 * Neo-Hookean laws are a subset)  */
#include "MooneyRivlinMaterialLaw.hpp"
/* Another incompressible material law */
#include "ExponentialMaterialLaw.hpp"
/* This is a useful helper class */
#include "NonlinearElasticityTools.hpp"
/* For visualising results in Paraview */
#include "VtkNonlinearElasticitySolutionWriter.hpp"

#include "PetscSetupUtils.hpp"
#include "CommandLineArguments.hpp"

#endif

namespace TimeSteppingScheme
{

template<int D>
QuasiStaticNonlinearElasticitySolverChaste<D>::
QuasiStaticNonlinearElasticitySolverChaste(DihuContext context) :
  context_(context["QuasiStaticNonlinearElasticitySolverChaste"]), data_(context_),
  initialized_(false)
{
  // get python config
  this->specificSettings_ = this->context_.getPythonConfig();

  // read in the durationLogKey
  if (specificSettings_.hasKey("durationLogKey"))
  {
    this->durationLogKey_ = specificSettings_.getOptionString("durationLogKey", "");
  }

  maximumActiveStress_ = specificSettings_.getOptionDouble("maximumActiveStress", 1.0, PythonUtility::ValidityCriterion::Positive);
  strainScalingCurveWidth_ = specificSettings_.getOptionDouble("strainScalingCurveWidth", 1.0, PythonUtility::ValidityCriterion::Positive);

  LOG(DEBUG) << "QuasiStaticNonlinearElasticitySolverChaste: parsed parameters maximumActiveStress: " << maximumActiveStress_ << ", strainScalingCurveWidth: " << strainScalingCurveWidth_;
  LOG(DEBUG) << "now parse output writers";

  // initialize output writers
  this->outputWriterManager_.initialize(this->context_, this->specificSettings_);
}

template<int D>
void QuasiStaticNonlinearElasticitySolverChaste<D>::
advanceTimeSpan()
{
  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  LOG(TRACE) << "advanceTimeSpan, endTime: " << endTime_;

  // compute active stress from activation

  LOG(DEBUG) << "solve linear elasticity";

#ifdef HAVE_CHASTE
  try
  {

    /* First, define the geometry. This should be specified using the `QuadraticMesh` class, which inherits from `TetrahedralMesh`
      * and has mostly the same interface. Here we define a 0.8 by 1 rectangle, with elements 0.1 wide.
      * (`QuadraticMesh`s can also be read in using `TrianglesMeshReader`; see next tutorial/rest of code base for examples of this).
      */
    //DistributedQuadraticMesh
    QuadraticMesh<3> mesh;
    global_no_t nElementsGlobalX = functionSpace_->meshPartition()->nElementsGlobal(0);
    global_no_t nElementsGlobalY = functionSpace_->meshPartition()->nElementsGlobal(1);
    global_no_t nElementsGlobalZ = functionSpace_->meshPartition()->nElementsGlobal(2);
    mesh.ConstructRegularSlabMesh(1.0 /*spaceStep*/, /*width*/ (double)nElementsGlobalX, /*height*/ (double)nElementsGlobalY, /*depth*/ (double)nElementsGlobalZ);

    LOG(DEBUG) << "chaste mesh has " << mesh.GetNumNodes() << " nodes, " << mesh.GetNumVertices() << " vertices, "
      << mesh.GetNumElements() << " elements.";

    // chaste mesh has 10309 nodes, 1519 vertices, 6480 elements.
    // note: mesh has tetrahedral elements, number of tretrahedral elements = 6*number of quadrilateral elements, node ordering?

    // find node indices where component 1 (y) has value
    std::vector<unsigned> topNodes = NonlinearElasticityTools<3>::GetNodesByComponentValue(mesh, /*component*/ 2, /*value*/ (double)nElementsGlobalZ);
    std::vector<unsigned> bottomNodes = NonlinearElasticityTools<3>::GetNodesByComponentValue(mesh, /*component*/ 2, /*value*/ 0.0);

    LOG(DEBUG) << "topNodes: " << topNodes;


    if (functionSpace_->meshPartition()->nNodesGlobal() != mesh.GetNumNodes())
    {
      LOG(ERROR) << "Number of nodes in opendihu parsed mesh (" << functionSpace_->meshPartition()->nNodesGlobal() << ") does not match "
        << "number of nodes in chaste mesh (" << mesh.GetNumNodes() << "). This means that transfer of the mesh from opendihu to chaste did not work.";
    }

    unsigned int nodeIndex = 0;
    c_vector<double,3> &location0 = mesh.GetNode(nodeIndex)->rGetModifiableLocation();
    location0[1] += 0.01;

    nodeIndex = 1;
    c_vector<double,3> &location1 = mesh.GetNode(nodeIndex)->rGetModifiableLocation();
    location1[1] += 0.02;

    nodeIndex = 2;
    c_vector<double,3> &location2 = mesh.GetNode(nodeIndex)->rGetModifiableLocation();
    location2[1] += 0.03;


    nodeIndex = 3;
    c_vector<double,3> &location3 = mesh.GetNode(nodeIndex)->rGetModifiableLocation();
    location3[1] += 0.11;

    nodeIndex = 4;
    c_vector<double,3> &location4 = mesh.GetNode(nodeIndex)->rGetModifiableLocation();
    location4[1] += 0.12;

    nodeIndex = 5;
    c_vector<double,3> &location5 = mesh.GetNode(nodeIndex)->rGetModifiableLocation();
    location5[1] += 0.13;

    // set new posititons
#if 1
    // numbering for chaste as shown in plot_chaste_numbering.py
    unsigned int chasteNodeIndex = 0;

    // loop over element corner nodes
    for (unsigned int nodeIndexZ = 0; nodeIndexZ < functionSpace_->meshPartition()->nNodesGlobal(2); nodeIndexZ += 2)
    {
      for (unsigned int nodeIndexY = 0; nodeIndexY < functionSpace_->meshPartition()->nNodesGlobal(1); nodeIndexY += 2)
      {
        for (unsigned int nodeIndexX = 0; nodeIndexX < functionSpace_->meshPartition()->nNodesGlobal(0); nodeIndexX += 2)
        {
          int opendihuNodeIndex = nodeIndexZ * functionSpace_->meshPartition()->nNodesGlobal(0) * functionSpace_->meshPartition()->nNodesGlobal(1)
            + nodeIndexY * functionSpace_->meshPartition()->nNodesGlobal(1) + nodeIndexX;

          c_vector<double,3> &location = mesh.GetNode(chasteNodeIndex)->rGetModifiableLocation();
          Vec3 position = functionSpace_->getGeometry(opendihuNodeIndex);
          location[0] = position[0];
          location[1] = position[1];
          location[2] = position[2];
          LOG(DEBUG) << "corner, nodeIndex " << opendihuNodeIndex << ", " << chasteNodeIndex << " pos " << location;

          chasteNodeIndex++;
        }
      }
    }

    // loop over rest of nodes
    for (unsigned int nodeIndexZ = 0; nodeIndexZ < functionSpace_->meshPartition()->nNodesGlobal(2); nodeIndexZ++)
    {
      for (unsigned int nodeIndexY = 0; nodeIndexY < functionSpace_->meshPartition()->nNodesGlobal(1); nodeIndexY++)
      {
        for (unsigned int nodeIndexX = 0; nodeIndexX < functionSpace_->meshPartition()->nNodesGlobal(0); nodeIndexX++)
        {
          // do not consider element corner nodes
          if (nodeIndexZ % 2 == 0 && nodeIndexY % 2 == 0 && nodeIndexX % 2 == 0)
            continue;

          int opendihuNodeIndex = nodeIndexZ * functionSpace_->meshPartition()->nNodesGlobal(0) * functionSpace_->meshPartition()->nNodesGlobal(1)
            + nodeIndexY * functionSpace_->meshPartition()->nNodesGlobal(1) + nodeIndexX;

          c_vector<double,3> &location = mesh.GetNode(chasteNodeIndex)->rGetModifiableLocation();
          Vec3 position = functionSpace_->getGeometry(opendihuNodeIndex);
          location[0] = position[0];
          location[1] = position[1];
          location[2] = position[2];
          LOG(DEBUG) << "nodeIndex " << opendihuNodeIndex << "," << chasteNodeIndex << " pos " << location;

          chasteNodeIndex++;
        }
      }
    }


#endif


    LOG(DEBUG) << "node positions set";

    // recalculate mesh properties, like Jacobian
    mesh.RefreshMesh();

    LOG(DEBUG) << "refreshed mesh";

    /* We use a Mooney-Rivlin material law, which applies to isotropic materials and has two parameters.
      * Restricted to 2D however, it only has one parameter, which can be thought of as the total
      * stiffness. We declare a Mooney-Rivlin law, setting the parameter to 1.
      */
    MooneyRivlinMaterialLaw<3> law(1.0, 0.5);

    /* Next, the body force density. In realistic problems this will either be
      * acceleration due to gravity (ie b=(0,-9.81)) or zero if the effect of gravity can be neglected.
      * In this problem we apply a gravity-like downward force.
      */
    c_vector<double,3> bodyForce;
    bodyForce(0) = 0.0;
    bodyForce(1) = 0.0;
    bodyForce(2) = -0.1;

    std::vector<BoundaryElement<2,3>*> tractionElements;
    std::vector<c_vector<double,3> > tractionValues;

    // Create a constant traction
    c_vector<double,3> traction;
    traction(0) = 0;
    traction(1) = 0;
    traction(2) = -0.1; // negative = outwards

    ChasteCuboid<3> boundingBox = mesh.CalculateBoundingBox();
    LOG(DEBUG) << "get bounding box: " << boundingBox.rGetLowerCorner().rGetLocation() << ", " << boundingBox.rGetUpperCorner().rGetLocation();

    int observationNodeNo = 0;

    // loop over boundary elements
    for (TetrahedralMesh<3,3>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
        iter != mesh.GetBoundaryElementIteratorEnd();
        ++iter)
    {
      LOG(DEBUG) << "boundary element centroid z: " << (*iter)->CalculateCentroid()[2];

      // if element lies at bottom of domain
      if (fabs((*iter)->CalculateCentroid()[2] - boundingBox.rGetLowerCorner().rGetLocation()[2]) < 1e-6)
      {
          // Put the boundary element and the constant traction into the stores.
          BoundaryElement<2,3>* p_element = *iter;
          tractionElements.push_back(p_element);
          tractionValues.push_back(traction);

          observationNodeNo = p_element->GetNodeGlobalIndex(0);
      }
    }



    /* Two types of boundary condition are required: displacement and traction. As with the other PDE solvers,
      * the displacement (Dirichlet) boundary conditions are specified at nodes, whereas traction (Neumann) boundary
      * conditions are specified on boundary elements.
      *
      * In this test we apply displacement boundary conditions on one surface of the mesh, the upper (Y=1.0) surface.
      * We are going to specify zero-displacement at these nodes.
      * We do not specify any traction boundary conditions, which means that (effectively) zero-traction boundary
      * conditions (ie zero pressures) are applied on the three other surfaces.
      *
      * We need to get a `std::vector` of all the node indices that we want to fix. The `NonlinearElasticityTools`
      * has a static method for helping do this: the following gets all the nodes for which Y=1.0. The second
      * argument (the '1') indicates Y . (So, for example, `GetNodesByComponentValue(mesh, 0, 10)` would get the nodes on X=10).
      */

    /*
      * Before creating the solver we create a `SolidMechanicsProblemDefinition` object,  which contains
      * everything that defines the problem: mesh, material law, body force,
      * the fixed nodes and their locations, any traction boundary conditions, and the density
      * (which multiplies the body force, otherwise isn't used).
      */

    SolidMechanicsProblemDefinition<3> problemDefinition(mesh);

    /* Set the material problem on the problem definition object, saying that the problem, and
      * the material law, is incompressible. All material law files can be found in
      * `continuum_mechanics/src/problem/material_laws`. */
    problemDefinition.SetMaterialLaw(INCOMPRESSIBLE,&law);

    /* Set the fixed nodes, choosing zero displacement for these nodes (see later for how
      * to provide locations for the fixed nodes). */
    problemDefinition.SetZeroDisplacementNodes(topNodes);
    /* Set the body force and the density. (Note that the second line isn't technically
      * needed, as internally the density is initialised to 1)
      */
    //problemDefinition.SetBodyForce(bodyForce);
    problemDefinition.SetDensity(1.0);
    problemDefinition.SetTractionBoundaryConditions(tractionElements, tractionValues);
    //problemDefinition.SetSolveUsingSnes(true);  // use PETSc SNES instead of own solver

    LOG(DEBUG) << tractionElements.size() << "traction elements: ";


    /* Now we create the (incompressible) solver, passing in the mesh, problem definition
      * and output directory, relative to $CHASTE_TEST_OUTPUT
      */
    IncompressibleNonlinearElasticitySolver<3> solver(mesh, problemDefinition, "output_directory");

    // enable active stress contributions
    solver.SetIncludeActiveTension(true);
    solver.SetTakeFullFirstNewtonStep(true);

    for (int loadStepIndex = 0; loadStepIndex < 10; loadStepIndex++)
    {

      solver.SetCurrentTime(loadStepIndex);

      // update traction
      c_vector<double,3> traction;
      traction(0) = 0;
      traction(1) = 0;
      traction(2) = -0.1*loadStepIndex; // negative = outwards

      // update the constant traction values
      tractionValues.clear();
      //std::fill_n(tractionValues.begin(), tractionElements.size(), traction);0
      for(int i = 0; i < tractionElements.size(); i++)
        tractionValues.push_back(traction);

      LOG(DEBUG) << "n traction values: " << tractionValues.size() << ", n traction elements: " << tractionElements.size();

      problemDefinition.SetTractionBoundaryConditions(tractionElements, tractionValues);

      LOG(DEBUG) << "solve...";

      /* .. and to compute the solution, just call `Solve()` */
      solver.Solve();
      LOG(DEBUG) << "solving done";

      /* '''Visualisation'''. Go to the folder `SimpleIncompressibleElasticityTutorial` in your test-output directory.
        * There should be 2 files, initial.nodes and solution.nodes. These are the original nodal positions and the deformed
        * positions. Each file has two columns, the x and y locations of each node. To visualise the solution in say
        * Matlab or Octave, you could do: `x=load('solution.nodes'); plot(x(:,1),x(:,2),'k*')`. For Cmgui output, see below.
        *
        * To get the actual solution from the solver, use these two methods. Note that the first
        * gets the deformed position (ie the new location, not the displacement). They are both of size
        * num_total_nodes.
        */
      std::vector<c_vector<double,3> >& r_deformed_positions = solver.rGetDeformedPosition();
      std::vector<double>& r_pressures = solver.rGetPressures();

      /* Let us obtain the values of the new position, and the pressure, at the bottom right corner node. */
      //assert( fabs(mesh.GetNode(node_index)->rGetLocation()[0] - 0.8) < 1e-6); // check that X=0.8, ie that we have the correct node,
      //assert( fabs(mesh.GetNode(node_index)->rGetLocation()[1] - 0.0) < 1e-6); // check that Y=0.0, ie that we have the correct node,
      std::cout << "New position: " << r_deformed_positions[observationNodeNo](0) << " " << r_deformed_positions[observationNodeNo](1) << " " << r_deformed_positions[observationNodeNo](2) << "\n";
      std::cout << "Pressure: " << r_pressures[observationNodeNo] << "\n";

      /* HOW_TO_TAG Continuum mechanics
        * Visualise nonlinear elasticity problems solutions, including visualisng strains
        */

      /* One visualiser is Cmgui. This method can be used to convert all the output files to Cmgui format.
        * They are placed in `[OUTPUT_DIRECTORY]/cmgui`. A script is created to easily load the data: in a
        * terminal cd to this directory and call `cmgui LoadSolutions.com`. (In this directory, the initial position is given by
        * solution_0.exnode, the deformed by solution_1.exnode).
        */
      //solver.CreateCmguiOutput();

      /* The recommended visualiser is Paraview, for which Chaste must be installed with VTK. With paraview, strains (and in the future
        * stresses) can be visualised on the undeformed/deformed geometry). We can create VTK output using
        * the `VtkNonlinearElasticitySolutionWriter` class. The undeformed geometry, solution displacement, and pressure (if incompressible
        * problem) are written to file, and below we also choose to write the deformation tensor C for each element.
        */
      VtkNonlinearElasticitySolutionWriter<3> vtk_writer(solver);
      vtk_writer.SetWriteElementWiseStrains(DEFORMATION_TENSOR_C); // other options are DEFORMATION_GRADIENT_F and LAGRANGE_STRAIN_E
      vtk_writer.Write();
    }
  }
  catch(Exception e)
  {
    LOG(ERROR) << "chaste exception: " << e.GetMessage();
  }
#endif


  // solve linear elasticity

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);

  // write current output values
  this->outputWriterManager_.writeOutput(this->data_, 0, endTime_);
}

template<int D>
void QuasiStaticNonlinearElasticitySolverChaste<D>::
run()
{
  // initialize everything
  initialize();

  this->advanceTimeSpan();
}

template<int D>
void QuasiStaticNonlinearElasticitySolverChaste<D>::
setTimeSpan(double startTime, double endTime)
{
  endTime_ = endTime;
}

template<int D>
void QuasiStaticNonlinearElasticitySolverChaste<D>::
initialize()
{
  if (this->initialized_)
    return;

  LOG(DEBUG) << "initialize QuasiStaticNonlinearElasticitySolverChaste";
  assert(this->specificSettings_.pyObject());

  // create function space / mesh, the geometry is from the settings
  LOG(DEBUG) << "FiniteElementMethodBase constructor, create new function space from settings";
  functionSpace_ = context_.meshManager()->functionSpace<FunctionSpaceType>(specificSettings_);

  // initialize the data object
  // store mesh in data
  data_.setFunctionSpace(functionSpace_);

  data_.initialize();

  LOG(DEBUG) << "initialization done";
  this->initialized_ = true;
}

template<int D>
void QuasiStaticNonlinearElasticitySolverChaste<D>::reset()
{
  this->initialized_ = false;
}

template<int D>
typename QuasiStaticNonlinearElasticitySolverChaste<D>::Data &QuasiStaticNonlinearElasticitySolverChaste<D>::
data()
{
  return data_;
}

//! get the data that will be transferred in the operator splitting to the other term of the splitting
//! the transfer is done by the slot_connector_data_transfer class
template<int D>
std::shared_ptr<typename QuasiStaticNonlinearElasticitySolverChaste<D>::SlotConnectorDataType>
QuasiStaticNonlinearElasticitySolverChaste<D>::
getSlotConnectorData()
{
  return nullptr;
}

//! output the given data for debugging
template<int D>
std::string QuasiStaticNonlinearElasticitySolverChaste<D>::
getString(std::shared_ptr<typename QuasiStaticNonlinearElasticitySolverChaste<D>::SlotConnectorDataType> data)
{
  std::stringstream s;
  s << "<QuasiStaticNonlinearElasticitySolverChaste:" << data << ">";
  return s.str();
}

} // namespace TimeSteppingScheme
