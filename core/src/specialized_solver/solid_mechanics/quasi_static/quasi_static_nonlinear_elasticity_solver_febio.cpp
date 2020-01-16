#include "specialized_solver/solid_mechanics/quasi_static/quasi_static_nonlinear_elasticity_solver_febio.h"

#include <Python.h>  // has to be the first included header

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "data_management/specialized_solver/multidomain.h"
#include "control/performance_measurement.h"

namespace TimeSteppingScheme
{

QuasiStaticNonlinearElasticitySolverFebio::
QuasiStaticNonlinearElasticitySolverFebio(DihuContext context) :
  context_(context["QuasiStaticNonlinearElasticitySolverFebio"]), data_(context_), initialized_(false)
{

  // get python config
  this->specificSettings_ = this->context_.getPythonConfig();

  // read in the durationLogKey
  if (specificSettings_.hasKey("durationLogKey"))
  {
    this->durationLogKey_ = specificSettings_.getOptionString("durationLogKey", "");
  }

  // load settings
  this->activationFactor_ = this->specificSettings_.getOptionDouble("activationFactor", 1e-5, PythonUtility::NonNegative);
  this->preLoadFactor_ = this->specificSettings_.getOptionDouble("preLoadFactor", 100, PythonUtility::NonNegative);

  // initialize output writers
  this->outputWriterManager_.initialize(this->context_, this->specificSettings_);

  LOG(DEBUG) << "initialized QuasiStaticNonlinearElasticitySolverFebio";
}


void QuasiStaticNonlinearElasticitySolverFebio::
advanceTimeSpan()
{
  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  createFebioInputFile();

  runFebio();

  loadFebioOutputFile();

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);

  // write current output values
  this->outputWriterManager_.writeOutput(this->data_, 0, endTime_);
}

void QuasiStaticNonlinearElasticitySolverFebio::
runFebio()
{
  // only run febio on rank 0
  int ownRankNo = DihuContext::ownRankNoCommWorld();

  if (ownRankNo == 0)
  {

    // run febio simulation
    int returnValue = system("rm -f febio_input.log febio_stress_output.txt febio_geometry_output.txt; febio2 -i febio_input.feb > /dev/null");

    if (returnValue != EXIT_SUCCESS)
    {
      LOG(ERROR) << "Running febio failed with error code " << returnValue;
    }

#ifndef NDEBUG
    // copy contents of log file from the febio execution to the opendihu.log file
    std::ifstream logFile("febio_input.log");

    if (logFile.is_open())
    {
      std::string logContents((std::istreambuf_iterator<char>(logFile)), std::istreambuf_iterator<char>());
      LOG(DEBUG) << "Content of febio log: " << std::endl << logContents;
    }
#endif
  }

  // barrier to wait until febio simulation is finished
  MPIUtility::handleReturnValue(MPI_Barrier(this->data_.functionSpace()->meshPartition()->mpiCommunicator()), "MPI_Barrier");
}

void QuasiStaticNonlinearElasticitySolverFebio::
createFebioInputFile()
{
  LOG(DEBUG) << "createFebioInputFile";

  int ownRankNo = DihuContext::ownRankNoCommWorld();

  // communicate elemental values
  std::vector<double> activationValuesGlobal;
  std::vector<int> nodeNosGlobal;
  communicateElementValues(activationValuesGlobal, nodeNosGlobal);

  // communicate node positions
  std::vector<double> nodePositionValuesGlobal;
  communicateNodeValues(nodePositionValuesGlobal);

  // loop over elements and create separate material for each element
  if (ownRankNo == 0)
  {
    std::stringstream fileContents;

    fileContents << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>" << std::endl
      << "<!--" << std::endl
      << "Problem Description:" << std::endl
      << "\tMuscle mesh to solve at time " << endTime_ << ", created by opendihu" << std::endl
      << "-->" << std::endl
      << "<febio_spec version=\"1.2\">" << std::endl
      << "\t<Control>" << std::endl
      << "\t\t<title>Biceps</title>" << std::endl
  //    << "\t\t<restart file=\"dump.out\">1</restart>" << std::endl
      << "\t\t<time_steps>10</time_steps>" << std::endl
      << "\t\t<step_size>0.1</step_size>" << std::endl
      << "\t\t<max_refs>15</max_refs>" << std::endl
      << "\t\t<max_ups>10</max_ups>" << std::endl
      << "\t\t<dtol>0.001</dtol>" << std::endl
      << "\t\t<etol>0.01</etol>" << std::endl
      << "\t\t<rtol>0</rtol>" << std::endl
      << "\t\t<lstol>0.9</lstol>" << std::endl
      << "\t\t<analysis type=\"dynamic\"></analysis>" << std::endl
      << "\t\t<time_stepper>" << std::endl
      << "\t\t\t<dtmin>0.01</dtmin>" << std::endl
      << "\t\t\t<dtmax>0.1</dtmax>" << std::endl
      << "\t\t\t<max_retries>5</max_retries>" << std::endl
      << "\t\t\t<opt_iter>10</opt_iter>" << std::endl
      << "\t\t</time_stepper>" << std::endl
      << "\t</Control>" << std::endl
      << "\t<Material>" << std::endl;

    // loop over global elements and add one material per element with the activation value of this element
    for (global_no_t elementNoGlobalPetsc = 0; elementNoGlobalPetsc < data_.functionSpace()->nElementsGlobal(); elementNoGlobalPetsc++)
    {
      double activationValue = activationValuesGlobal[elementNoGlobalPetsc];

      fileContents << "\t\t<material id=\"" << elementNoGlobalPetsc+1 << "\" name=\"Material for element " << elementNoGlobalPetsc << "\" type=\"muscle material\">" << std::endl
        << "\t\t\t<g1>500</g1>" << std::endl
        << "\t\t\t<g2>500</g2>" << std::endl
        << "\t\t\t<p1>0.05</p1>" << std::endl
        << "\t\t\t<p2>6.6</p2>" << std::endl
        << "\t\t\t<smax>3e5</smax>" << std::endl
        << "\t\t\t<Lofl>1.07</Lofl>" << std::endl
        << "\t\t\t<lam_max>1.4</lam_max>" << std::endl
        << "\t\t\t<k>1e6</k>" << std::endl
        << "\t\t\t<fiber type=\"local\">1, 5</fiber>    <!-- fiber direction is from local node 1 to 5 in every element, i.e. upwards -->" << std::endl
        << "\t\t\t<activation lc=\"1\">" << activationValue << "</activation>" << std::endl
        << "\t\t</material>" << std::endl;
    }

    fileContents << "\t</Material>" << std::endl
      << "\t<Geometry>" << std::endl
      << "\t\t<Nodes>" << std::endl;

    // loop over nodes
    for (global_no_t nodeNoGlobalPetsc = 0; nodeNoGlobalPetsc < data_.functionSpace()->nNodesGlobal(); nodeNoGlobalPetsc++)
    {
      double nodePositionX = nodePositionValuesGlobal[3*nodeNoGlobalPetsc + 0];
      double nodePositionY = nodePositionValuesGlobal[3*nodeNoGlobalPetsc + 1];
      double nodePositionZ = nodePositionValuesGlobal[3*nodeNoGlobalPetsc + 2];

      fileContents << "\t\t\t<node id=\"" << nodeNoGlobalPetsc+1 << "\">"
        << nodePositionX << "," << nodePositionY << "," << nodePositionZ << "</node>" << std::endl;
    }

    fileContents << "\t\t</Nodes>" << std::endl
      << "\t\t<Elements>" << std::endl;

    // loop over global elements and insert node nos of elements
    for (global_no_t elementNoGlobalPetsc = 0; elementNoGlobalPetsc < data_.functionSpace()->nElementsGlobal(); elementNoGlobalPetsc++)
    {
      nodeNosGlobal[elementNoGlobalPetsc*8 + 0];

      // elements types: https://help.febio.org/FEBio/FEBio_um_2_8/FEBio_um_2-8-3.8.2.1.html
      fileContents << "\t\t\t<hex8 id=\"" << elementNoGlobalPetsc+1 << "\" mat=\"" << elementNoGlobalPetsc+1 << "\"> ";
      fileContents << nodeNosGlobal[elementNoGlobalPetsc*8 + 0] << ", " << nodeNosGlobal[elementNoGlobalPetsc*8 + 1]
        << ", " << nodeNosGlobal[elementNoGlobalPetsc*8 + 2] << ", " << nodeNosGlobal[elementNoGlobalPetsc*8 + 3]
        << ", " << nodeNosGlobal[elementNoGlobalPetsc*8 + 4] << ", " << nodeNosGlobal[elementNoGlobalPetsc*8 + 5]
        << ", " << nodeNosGlobal[elementNoGlobalPetsc*8 + 6] << ", " << nodeNosGlobal[elementNoGlobalPetsc*8 + 7];
      fileContents << "</hex8>" << std::endl;
    }

    fileContents << "\t\t</Elements>" << std::endl
      << "\t</Geometry>" << std::endl
      << "\t<Boundary>" << std::endl
      << "\t\t<fix>" << std::endl;

    // fix bottom dofs
    int nNodesX = this->data_.functionSpace()->nNodesGlobal(0);
    int nNodesY = this->data_.functionSpace()->nNodesGlobal(1);
    int nNodesZ = this->data_.functionSpace()->nNodesGlobal(2);

    // fix x direction for left row
    for (int j = 0; j < nNodesY; j++)
    {
      dof_no_t dofNoLocal = j*nNodesX;
      fileContents << "\t\t\t<node id=\"" << dofNoLocal+1 << "\" bc=\"x\"></node>" << std::endl;
    }

    fileContents << "\t\t</fix>" << std::endl
      << "\t\t<fix> " << std::endl;

    // fix y direction for front row
    for (int i = 0; i < nNodesX; i++)
    {
      dof_no_t dofNoLocal = i;
      fileContents << "\t\t\t<node id=\"" << dofNoLocal+1 << "\" bc=\"y\"></node>" << std::endl;
    }


    fileContents << "\t\t</fix>" << std::endl
      << "\t\t<fix> " << std::endl;

    // fix z direction for all
    for (int j = 0; j < nNodesY; j++)
    {
      for (int i = 0; i < nNodesX; i++)
      {
        dof_no_t dofNoLocal = j*nNodesX + i;
        fileContents << "\t\t\t<node id=\"" << dofNoLocal+1 << "\" bc=\"z\"></node>" << std::endl;
      }
    }

    fileContents << "\t\t</fix>" << std::endl
      << "\t</Boundary>" << std::endl;

    // force that stretches the muscle
    fileContents << "\t<Loads>" << std::endl
      << "\t\t<force>" << std::endl;

    // add force for top nodes
    for (int j = 0; j < nNodesY; j++)
    {
      for (int i = 0; i < nNodesX; i++)
      {
        dof_no_t dofNoLocal = (nNodesZ-1)*nNodesY*nNodesX + j*nNodesX + i;

        std::string value = "1.0";
        if (i == 0 || i == nNodesX-1)
        {
          if (j == 0 || j == nNodesY-1)
          {
            value = "0.25";
          }
          else
          {
            value = "0.5";
          }
        }
        else if (j == 0 || j == nNodesY-1)
        {
          value = "0.5";
        }

        fileContents << "\t\t\t<node id=\"" << dofNoLocal+1 << "\" bc=\"z\" lc=\"2\">"
          << value << "</node>" << std::endl;
      }
    }


    fileContents << "\t\t</force>" << std::endl
      << "\t</Loads>" << std::endl;

    fileContents << "\t<LoadData>" << std::endl
      << "\t\t<loadcurve id=\"1\"> <!-- for activation -->" << std::endl
      << "\t\t\t<loadpoint>0,0</loadpoint>" << std::endl
      << "\t\t\t<loadpoint>1," << activationFactor_ << "</loadpoint>" << std::endl
      << "\t\t</loadcurve>" << std::endl
      << "\t\t<loadcurve id=\"2\"> <!-- for external load that stretches the muscle -->" << std::endl
      << "\t\t\t<loadpoint>0,0</loadpoint>" << std::endl
      << "\t\t\t<loadpoint>1," << preLoadFactor_ << "</loadpoint>" << std::endl
      << "\t\t</loadcurve>" << std::endl
      << "\t</LoadData>" << std::endl
      << "\t<Output>" << std::endl
      << "\t\t<plotfile type=\"febio\">" << std::endl
      << "\t\t\t<var type=\"fiber vector\"/>" << std::endl
      << "\t\t\t<var type=\"displacement\"/>" << std::endl
      << "\t\t\t<var type=\"fiber stretch\"/>" << std::endl
      << "\t\t\t<var type=\"Lagrange strain\"/>" << std::endl
      << "\t\t\t<var type=\"parameter\"/>" << std::endl
      << "\t\t\t<var type=\"relative volume\"/>" << std::endl
      << "\t\t\t<var type=\"stress\"/>" << std::endl
      << "\t\t\t<var type=\"strain energy density\"/>" << std::endl
      << "\t\t</plotfile>" << std::endl
      << "\t\t<logfile>" << std::endl

      // available variables: https://help.febio.org/FEBio/FEBio_um_2_9/index.html Sec. 3.17.1.2 and 3.17.1.3
      << "\t\t\t<node_data file=\"febio_geometry_output.txt\" format=\"%i,%g,%g,%g,%g,%g,%g,%g,%g,%g\" data=\"x;y;z;ux;uy;uz;Rx;Ry;Rz\"/>" << std::endl
      << "\t\t\t<element_data file=\"febio_stress_output.txt\" format=\"%i,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\" data=\"sx;sy;sz;sxy;syz;sxz;Ex;Ey;Ez;Exy;Eyz;Exz;J\"/>" << std::endl
      << "\t\t</logfile>" << std::endl
      << "\t</Output>" << std::endl
      << "</febio_spec>" << std::endl;


    std::ofstream file("febio_input.feb");

    if (!file.is_open())
    {
      LOG(ERROR) << "Could not write to file \"febio_input.feb\".";
    }

    file << fileContents.str();
    file.close();
  }
}

void QuasiStaticNonlinearElasticitySolverFebio::
loadFebioOutputFile()
{
  LOG(TRACE) << "loadFebioOutputFile";

  std::ifstream fileGeometry;
  fileGeometry.open("febio_geometry_output.txt", std::ios::in | std::ios::binary);

  if (!fileGeometry.is_open())
  {
    LOG(WARNING) << "Could not read output file \"febio_geometry_output.txt\" that should have been created by FEBio.";
    return;
  }

  // example of a file
  /*
*Step  = 11
*Time  = 1
*Data  = x;y;z
1,0.00284281,0,0
2,0.5,0,0
3,1,0,0
4,0,0.5,0
*/

  // determine number of time steps that are contained in file
  std::string line;
  int nStepsContainedInFile = 0;
  while(!fileGeometry.eof())
  {
    std::getline(fileGeometry, line);
    if (line.find("*Step") != std::string::npos)
    {
      nStepsContainedInFile = (int)(atoi(line.substr(line.find("=")+1).c_str()));
    }
  }

  LOG(DEBUG) << "nStepsContainedInFile: " << nStepsContainedInFile;

  // rewind file pointer
  fileGeometry.clear();
  fileGeometry.seekg(0);

  // load last step
  int currentStepNo = 0;
  std::vector<Vec3> geometryValues;
  while(!fileGeometry.eof())
  {
    std::getline(fileGeometry, line);
    if (line.find("*Step") != std::string::npos)
    {
      currentStepNo = (int)(atoi(line.substr(line.find("=")+1).c_str()));
    }

    if (currentStepNo == nStepsContainedInFile)
    {
      if (line.find("*") == std::string::npos && line.find(",") != std::string::npos)
      {
        int id = atoi(StringUtility::extractUntil(line, ",").c_str());

        // x;y;z;ux;uy;uz;Rx;Ry;Rz
        double x = atof(StringUtility::extractUntil(line, ",").c_str());
        double y = atof(StringUtility::extractUntil(line, ",").c_str());
        double z = atof(StringUtility::extractUntil(line, ",").c_str());
        double ux = atof(StringUtility::extractUntil(line, ",").c_str());
        double uy = atof(StringUtility::extractUntil(line, ",").c_str());
        double uz = atof(StringUtility::extractUntil(line, ",").c_str());
        double Rx = atof(StringUtility::extractUntil(line, ",").c_str());
        double Ry = atof(StringUtility::extractUntil(line, ",").c_str());
        double Rz = atof(line.c_str());

        global_no_t nodeNoGlobalPetsc = id - 1;

        // convert global node no to local no
        bool isLocal = false;
        node_no_t nodeNoLocal = this->data_.functionSpace()->meshPartition()->getNodeNoLocal(nodeNoGlobalPetsc, isLocal);

        LOG(DEBUG) << "read point " << id << ": " << Vec3({x,y,z}) << ", global: " << nodeNoGlobalPetsc << ", local: " << nodeNoLocal;

        if (isLocal)
        {
          LOG(DEBUG) << "is local";

          this->data_.displacements()->setValue(nodeNoLocal, Vec3{ux,uy,uz});
          this->data_.reactionForce()->setValue(nodeNoLocal, Vec3{Rx,Ry,Rz});

          geometryValues.push_back(Vec3{x,y,z});
        }
      }
    }
  }

  fileGeometry.close();

  std::ifstream fileStress;
  fileStress.open("febio_stress_output.txt", std::ios::in | std::ios::binary);

  if (!fileStress.is_open())
  {
    LOG(WARNING) << "Could not read output file \"febio_stress_output.txt\" that should have been created by FEBio.";
    return;
  }

  // load last step
  currentStepNo = 0;
  std::vector<Vec3> elementValues;

  // advance to line in file where actual data starts
  while(!fileStress.eof())
  {
    std::getline(fileStress, line);
    if (line.find("*Step") != std::string::npos)
    {
      currentStepNo = (int)(atoi(line.substr(line.find("=")+1).c_str()));
      LOG(DEBUG) << "currentStepNo: " << currentStepNo;
    }

    if (currentStepNo == nStepsContainedInFile)
    {
      if (line.find("*") == std::string::npos)
      {
        LOG(DEBUG) << "line to start: [" << line << "]";
        break;
      }
    }
  }

  this->data_.cauchyStress()->zeroEntries();
  this->data_.greenLagrangeStrain()->zeroEntries();
  this->data_.relativeVolume()->zeroEntries();

  std::vector<int> nSummands(this->data_.functionSpace()->nNodesLocalWithoutGhosts(), 0);
  int nElementsLoaded = 0;

  // loop over elements, average element-based values to nodal values
  for (;;)
  {
    int id = atoi(StringUtility::extractUntil(line, ",").c_str());

    // get local element no for global element no
    global_no_t elementNoGlobalPetsc = id - 1;
    bool isOnLocalDomain = false;
    element_no_t elementNoLocal = data_.functionSpace()->meshPartition()->getElementNoLocal(elementNoGlobalPetsc, isOnLocalDomain);

    LOG(DEBUG) << "read element global " << elementNoGlobalPetsc << ", local: " << elementNoLocal << ", isOnLocalDomain: " << isOnLocalDomain;

    if (isOnLocalDomain)
    {
      std::array<dof_no_t,FunctionSpace::nNodesPerElement()> elementNodeNos = data_.functionSpace()->getElementNodeNos(elementNoLocal);

      // sx;sy;sz;sxy;syz;sxz;Ex;Ey;Ez;Exy;Eyz;Exz;J
      double sx  = atof(StringUtility::extractUntil(line, ",").c_str());
      double sy  = atof(StringUtility::extractUntil(line, ",").c_str());
      double sz  = atof(StringUtility::extractUntil(line, ",").c_str());
      double sxy = atof(StringUtility::extractUntil(line, ",").c_str());
      double syz = atof(StringUtility::extractUntil(line, ",").c_str());
      double sxz = atof(StringUtility::extractUntil(line, ",").c_str());
      double Ex  = atof(StringUtility::extractUntil(line, ",").c_str());
      double Ey  = atof(StringUtility::extractUntil(line, ",").c_str());
      double Ez  = atof(StringUtility::extractUntil(line, ",").c_str());
      double Exy = atof(StringUtility::extractUntil(line, ",").c_str());
      double Eyz = atof(StringUtility::extractUntil(line, ",").c_str());
      double Exz = atof(StringUtility::extractUntil(line, ",").c_str());
      double J   = atof(line.c_str());

      LOG(DEBUG) << "local element " << elementNoLocal << " of " << this->data_.functionSpace()->nElementsLocal() << ", " << FunctionSpace::nNodesPerElement() << " elementNodeNos: " << elementNodeNos;

      for (node_no_t elementalNodeNo = 0; elementalNodeNo < FunctionSpace::nNodesPerElement(); elementalNodeNo++)
      {
        node_no_t nodeNoLocal = elementNodeNos[elementalNodeNo];

        if (nodeNoLocal >= this->data_.functionSpace()->nNodesLocalWithoutGhosts())
          continue;

        this->data_.cauchyStress()->setValue(nodeNoLocal, std::array<double,6>{sx,sy,sz,sxy,syz,sxz}, ADD_VALUES);
        this->data_.greenLagrangeStrain()->setValue(nodeNoLocal, std::array<double,6>{Ex,Ey,Ez,Exy,Eyz,Exz}, ADD_VALUES);
        this->data_.relativeVolume()->setValue(nodeNoLocal, J, ADD_VALUES);

        assert (nodeNoLocal < nSummands.size());
        nSummands[nodeNoLocal]++;
      }
      nElementsLoaded++;
    }

    // read next line from file
    std::getline(fileStress, line);

    // check if file is done
    if (line.find("*") != std::string::npos || fileStress.eof())
    {
      if (nElementsLoaded != data_.functionSpace()->nElementsLocal())
      {
        LOG(ERROR) << "Loaded " << nElementsLoaded << " local elements from \"febio_stress_output.txt\", but " << this->data_.functionSpace()->nElementsLocal() << " were expected.";
      }
      break;
    }
  }

  // divide values at nodes by number of summands
  // loop over nodes
  for (node_no_t nodeNoLocal = 0; nodeNoLocal < data_.functionSpace()->nNodesLocalWithoutGhosts(); nodeNoLocal++)
  {
    // cauchy stress
    std::array<double,6> cauchyStressValues = this->data_.cauchyStress()->getValue(nodeNoLocal);
    cauchyStressValues /= nSummands[nodeNoLocal];
    this->data_.cauchyStress()->setValue(nodeNoLocal, cauchyStressValues, INSERT_VALUES);

    // green lagrange strain
    std::array<double,6> greenLagrangeStrainValues = this->data_.greenLagrangeStrain()->getValue(nodeNoLocal);
    greenLagrangeStrainValues /= nSummands[nodeNoLocal];
    this->data_.greenLagrangeStrain()->setValue(nodeNoLocal, greenLagrangeStrainValues, INSERT_VALUES);

    // relative volume
    double relativeVolumeValue = this->data_.relativeVolume()->getValue(nodeNoLocal);
    relativeVolumeValue /= nSummands[nodeNoLocal];
    this->data_.relativeVolume()->setValue(nodeNoLocal, relativeVolumeValue, INSERT_VALUES);
  }

  fileStress.close();

  // update function space
  LOG(DEBUG) << "geometry field has representation "
    << this->data_.functionSpace()->geometryField().partitionedPetscVec()->getCurrentRepresentationString();

  this->data_.functionSpace()->geometryField().finishGhostManipulation();
  this->data_.functionSpace()->geometryField().setValuesWithoutGhosts(geometryValues);

  this->data_.functionSpace()->geometryField().zeroGhostBuffer();
  this->data_.functionSpace()->geometryField().finishGhostManipulation();
  this->data_.functionSpace()->geometryField().startGhostManipulation();

  LOG(DEBUG) << "geometryField pointer: " << this->data_.functionSpace()->geometryField().partitionedPetscVec();
  LOG(DEBUG) << "referenceGeometry pointer: " << this->data_.referenceGeometry()->partitionedPetscVec();
}

void QuasiStaticNonlinearElasticitySolverFebio::
communicateElementValues(std::vector<double> &activationValuesGlobal, std::vector<int> &nodeNosGlobal)
{
  int ownRankNo = DihuContext::ownRankNoCommWorld();
  int nRanks = DihuContext::nRanksCommWorld();

  // collect local activation values and node nos
  int ownSize = data_.functionSpace()->nElementsLocal();
  std::vector<double> activationValuesLocal(ownSize);
  std::vector<int> nodeNosLocal(8*ownSize);

  for (element_no_t elementNoLocal = 0; elementNoLocal < ownSize; elementNoLocal++)
  {
    std::array<double,FunctionSpace::nDofsPerElement()> activationElementalValues;
    data_.activation()->getElementValues(elementNoLocal, activationElementalValues);

    LOG(DEBUG) << "element " << elementNoLocal << ", activationElementalValues: " << activationElementalValues;

    Vec3 xi({0.5,0.5,0.5});
    double activationValue = this->data().functionSpace()->interpolateValueInElement(activationElementalValues, xi);

    activationValuesLocal[elementNoLocal] = activationValue;


    std::array<dof_no_t,FunctionSpace::nNodesPerElement()> elementNodeNos = data_.functionSpace()->getElementNodeNos(elementNoLocal);

    std::vector<dof_no_t> dofNosLocal(elementNodeNos.begin(), elementNodeNos.end());
    std::vector<PetscInt> dofNosGlobalPetsc;
    this->data().functionSpace()->meshPartition()->getDofNoGlobalPetsc(dofNosLocal, dofNosGlobalPetsc);

    nodeNosLocal[elementNoLocal*8 + 0] = dofNosGlobalPetsc[0]+1;
    nodeNosLocal[elementNoLocal*8 + 1] = dofNosGlobalPetsc[1]+1;
    nodeNosLocal[elementNoLocal*8 + 2] = dofNosGlobalPetsc[3]+1;
    nodeNosLocal[elementNoLocal*8 + 3] = dofNosGlobalPetsc[2]+1;
    nodeNosLocal[elementNoLocal*8 + 4] = dofNosGlobalPetsc[4]+1;
    nodeNosLocal[elementNoLocal*8 + 5] = dofNosGlobalPetsc[5]+1;
    nodeNosLocal[elementNoLocal*8 + 6] = dofNosGlobalPetsc[7]+1;
    nodeNosLocal[elementNoLocal*8 + 7] = dofNosGlobalPetsc[6]+1;
  }

  // perpare helper values for MPI_Gatherv
  std::vector<int> sizesOnRanks(nRanks);

  //std::cout << "at StimulationLogging::logStimulationBegin, ownRankNo: " << ownRankNo << ", nRanks: " << nRanks << std::endl;
  MPIUtility::handleReturnValue(MPI_Allgather(&ownSize, 1, MPI_INT, sizesOnRanks.data(), 1, MPI_INT, MPI_COMM_WORLD), "MPI_Allgather");

  std::vector<int> sizesOnRanksNodeNos(nRanks);
  for (std::vector<int>::const_iterator iter = sizesOnRanks.cbegin(); iter != sizesOnRanks.cend(); iter++)
  {
    sizesOnRanksNodeNos[iter - sizesOnRanks.begin()] = (*iter) * 8;
  }

  // count total number of entries
  int nTotalEntries = 0;
  for (int rankNo = 0; rankNo < nRanks; rankNo++)
  {
    nTotalEntries += sizesOnRanks[rankNo];
  }

  assert(nTotalEntries == data_.functionSpace()->nElementsGlobal());

  LOG(DEBUG) << "nTotalEntries: " << nTotalEntries << ", sizesOnRanks: " << sizesOnRanks;

  if (nTotalEntries == 0)
    return;

  // communicate all values to rank 0
  // setup offsets for MPI_Gatherv
  std::vector<int> offsets(nRanks);
  std::vector<int> offsetsNodeNos(nRanks);

  for (int rankNo = 0; rankNo < nRanks; rankNo++)
  {
    if (rankNo == 0)
    {
      offsets[rankNo] = 0;
      offsetsNodeNos[rankNo] = 0;
    }
    else
    {
      offsets[rankNo] = offsets[rankNo-1] + sizesOnRanks[rankNo-1];
      offsetsNodeNos[rankNo] = 8*offsets[rankNo];
    }
  }

  // communicate activation values
  if (ownRankNo == 0)
  {
    activationValuesGlobal.resize(nTotalEntries);
    nodeNosGlobal.resize(8*nTotalEntries);
  }

  MPIUtility::handleReturnValue(MPI_Gatherv(activationValuesLocal.data(), sizesOnRanks[ownRankNo], MPI_DOUBLE,
                                            activationValuesGlobal.data(), sizesOnRanks.data(), offsets.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD), "MPI_Gatherv");
  MPIUtility::handleReturnValue(MPI_Gatherv(nodeNosLocal.data(),  sizesOnRanksNodeNos[ownRankNo], MPI_INT,
                                            nodeNosGlobal.data(), sizesOnRanksNodeNos.data(), offsetsNodeNos.data(), MPI_INT, 0, MPI_COMM_WORLD), "MPI_Gatherv");

  if (ownRankNo == 0)
  {
    LOG(DEBUG) << "activationValuesLocal: " << activationValuesLocal << ", activationValuesGlobal: " << activationValuesGlobal;
  }
}

void QuasiStaticNonlinearElasticitySolverFebio::
communicateNodeValues(std::vector<double> &nodePositionValuesGlobal)
{
  int ownRankNo = DihuContext::ownRankNoCommWorld();
  int nRanks = DihuContext::nRanksCommWorld();

  // collect local node positions
  int ownSize = data_.functionSpace()->nNodesLocalWithoutGhosts() * 3;
  std::vector<double> nodePositionValuesLocal(ownSize);

  // loop over nodes
  for (node_no_t nodeNoLocal = 0; nodeNoLocal < data_.functionSpace()->nNodesLocalWithoutGhosts(); nodeNoLocal++)
  {
    dof_no_t dofNoLocal = nodeNoLocal;
    Vec3 position = data_.referenceGeometry()->getValue(dofNoLocal);

    for (int i = 0; i < 3; i++)
      nodePositionValuesLocal[3*dofNoLocal + i] = position[i];
  }

  // perpare helper values for MPI_Gatherv
  std::vector<int> sizesOnRanks(nRanks);

  //std::cout << "at StimulationLogging::logStimulationBegin, ownRankNo: " << ownRankNo << ", nRanks: " << nRanks << std::endl;
  MPIUtility::handleReturnValue(MPI_Allgather(&ownSize, 1, MPI_INT, sizesOnRanks.data(), 1, MPI_INT, MPI_COMM_WORLD), "MPI_Allgather");

  // count total number of entries
  int nTotalEntries = 0;
  for (int rankNo = 0; rankNo < nRanks; rankNo++)
  {
    nTotalEntries += sizesOnRanks[rankNo];
  }

  assert(nTotalEntries == 3*data_.functionSpace()->nNodesGlobal());

  LOG(DEBUG) << "nTotalEntries: " << nTotalEntries << ", sizesOnRanks: " << sizesOnRanks;

  if (nTotalEntries == 0)
    return;

  // communicate all values to rank 0
  // setup offsets for MPI_Gatherv
  std::vector<int> offsets(nRanks);

  for (int rankNo = 0; rankNo < nRanks; rankNo++)
  {
    if (rankNo == 0)
    {
      offsets[rankNo] = 0;
    }
    else
    {
      offsets[rankNo] = offsets[rankNo-1] + sizesOnRanks[rankNo-1];
    }
  }

  // communicate activation values
  if (ownRankNo == 0)
    nodePositionValuesGlobal.resize(nTotalEntries);

  MPIUtility::handleReturnValue(MPI_Gatherv(nodePositionValuesLocal.data(), sizesOnRanks[ownRankNo], MPI_DOUBLE,
                                            nodePositionValuesGlobal.data(), sizesOnRanks.data(), offsets.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD), "MPI_Gatherv");
}

void QuasiStaticNonlinearElasticitySolverFebio::
run()
{
  // initialize everything
  initialize();

  this->advanceTimeSpan();
}


void QuasiStaticNonlinearElasticitySolverFebio::
setTimeSpan(double startTime, double endTime)
{
  endTime_ = endTime;
}


void QuasiStaticNonlinearElasticitySolverFebio::
initialize()
{
  LOG(DEBUG) << "initialize QuasiStaticNonlinearElasticitySolverFebio";
  assert(this->specificSettings_.pyObject());

  // create function space
  std::shared_ptr<FunctionSpace> functionSpace = context_.meshManager()->functionSpace<FunctionSpace>(specificSettings_);

  // initialize the data object
  // store mesh in data
  data_.setFunctionSpace(functionSpace);
  data_.initialize();

  // write initial geometry
  this->outputWriterManager_.writeOutput(this->data_, 0, 0);

  LOG(DEBUG) << "initialization done";
  this->initialized_ = true;
}


void QuasiStaticNonlinearElasticitySolverFebio::reset()
{
  this->initialized_ = false;
}


typename QuasiStaticNonlinearElasticitySolverFebio::Data &QuasiStaticNonlinearElasticitySolverFebio::
data()
{
  return data_;
}

//! get the data that will be transferred in the operator splitting to the other term of the splitting
//! the transfer is done by the output_connector_data_transfer class

std::shared_ptr<typename QuasiStaticNonlinearElasticitySolverFebio::OutputConnectorDataType>
QuasiStaticNonlinearElasticitySolverFebio::
getOutputConnectorData()
{
  return this->data_.getOutputConnectorData();
}

//! output the given data for debugging

std::string QuasiStaticNonlinearElasticitySolverFebio::
getString(std::shared_ptr<typename QuasiStaticNonlinearElasticitySolverFebio::OutputConnectorDataType> data)
{
  std::stringstream s;
  //s << "<QuasiStaticNonlinearElasticitySolverFebio:" << *data.activation() << ">";
  return s.str();
}

} // namespace TimeSteppingScheme
