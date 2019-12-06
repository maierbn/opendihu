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
  // run febio
  int returnValue = system("rm -f febio_input.log febio_stress_output.txt febio_geometry_output.txt; febio2 -i febio_input.feb > /dev/null");

  if (returnValue != EXIT_SUCCESS)
  {
    LOG(ERROR) << "Running febio failed with error code " << returnValue;
  }

#ifndef NDEBUG
  std::ifstream logFile("febio_input.log");

  if (logFile.is_open())
  {
    std::string logContents((std::istreambuf_iterator<char>(logFile)), std::istreambuf_iterator<char>());
    LOG(DEBUG) << "Content of febio log: " << std::endl << logContents;
  }
#endif
}

void QuasiStaticNonlinearElasticitySolverFebio::
createFebioInputFile()
{
  LOG(DEBUG) << "createFebioInputFile";

  std::ofstream file("febio_input.feb");

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

  // loop over elements and create separate material for each element
  for (element_no_t elementNoLocal = 0; elementNoLocal < data_.functionSpace()->nElementsLocal(); elementNoLocal++)
  {
    std::array<double,FunctionSpace::nDofsPerElement()> activationElementalValues;
    data_.activation()->getElementValues(elementNoLocal, activationElementalValues);
    Vec3 xi({0.5,0.5,0.5});
    double activationValue = this->data().functionSpace()->interpolateValueInElement(activationElementalValues, xi);

    fileContents << "\t\t<material id=\"" << elementNoLocal+1 << "\" name=\"Material for elmement " << elementNoLocal << "\" type=\"muscle material\">" << std::endl
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
  for (node_no_t nodeNoLocal = 0; nodeNoLocal < data_.functionSpace()->nNodesLocalWithoutGhosts(); nodeNoLocal++)
  {
    dof_no_t dofNoLocal = nodeNoLocal;
    //Vec3 position = data_.functionSpace()->geometryField().getValue(dofNoLocal);
    Vec3 position = data_.referenceGeometry()->getValue(dofNoLocal);
    fileContents << "\t\t\t<node id=\"" << nodeNoLocal+1 << "\">" << position[0] << "," << position[1] << "," << position[2]
      << "</node>" << std::endl;
  }

  fileContents << "\t\t</Nodes>" << std::endl
    << "\t\t<Elements>" << std::endl;

  // loop over elements
  for (element_no_t elementNoLocal = 0; elementNoLocal < data_.functionSpace()->nElementsLocal(); elementNoLocal++)
  {
    std::array<dof_no_t,FunctionSpace::nNodesPerElement()> elementNodeNos = data_.functionSpace()->getElementNodeNos(elementNoLocal);

    // elements types: https://help.febio.org/FEBio/FEBio_um_2_8/FEBio_um_2-8-3.8.2.1.html
    fileContents << "\t\t\t<hex8 id=\"" << elementNoLocal+1 << "\" mat=\"" << elementNoLocal+1 << "\"> ";
    fileContents << elementNodeNos[0]+1 << ", " << elementNodeNos[1]+1 << ", " << elementNodeNos[3]+1 << ", " << elementNodeNos[2]+1 << ", "
      << elementNodeNos[4]+1 << ", " << elementNodeNos[5]+1 << ", " << elementNodeNos[7]+1 << ", " << elementNodeNos[6]+1;
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
    << "\t\t\t<loadpoint>1,1e-5</loadpoint>" << std::endl
    << "\t\t</loadcurve>" << std::endl
		<< "\t\t<loadcurve id=\"2\"> <!-- for external load that stretches the muscle -->" << std::endl
    << "\t\t\t<loadpoint>0,0</loadpoint>" << std::endl
		<< "\t\t\t<loadpoint>1,100</loadpoint>" << std::endl
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

  file << fileContents.str();
  file.close();
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

        LOG(DEBUG) << "read point " << id << ": " << Vec3({x,y,z});
        geometryValues.push_back(Vec3{x,y,z});

        node_no_t nodeNo = id - 1;
        this->data_.displacements()->setValue(nodeNo, Vec3{ux,uy,uz});
        this->data_.reactionForce()->setValue(nodeNo, Vec3{Rx,Ry,Rz});
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

  // loop over elements, average element-based values to nodal values
  for (element_no_t elementNoLocal = 0; elementNoLocal < data_.functionSpace()->nElementsLocal(); elementNoLocal++)
  {
    std::array<dof_no_t,FunctionSpace::nNodesPerElement()> elementNodeNos = data_.functionSpace()->getElementNodeNos(elementNoLocal);

    int id = atoi(StringUtility::extractUntil(line, ",").c_str());
    LOG(DEBUG) << "id: " << id;

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

    for (node_no_t elementalNodeNo = 0; elementalNodeNo < FunctionSpace::nNodesPerElement(); elementalNodeNo++)
    {
      node_no_t nodeNoLocal = elementNodeNos[elementalNodeNo];
      this->data_.cauchyStress()->setValue(nodeNoLocal, std::array<double,6>{sx,sy,sz,sxy,syz,sxz}, ADD_VALUES);
      this->data_.greenLagrangeStrain()->setValue(nodeNoLocal, std::array<double,6>{Ex,Ey,Ez,Exy,Eyz,Exz}, ADD_VALUES);
      this->data_.relativeVolume()->setValue(nodeNoLocal, J, ADD_VALUES);
      nSummands[nodeNoLocal]++;
    }

    std::getline(fileStress, line);
    if (line.find("*") != std::string::npos || fileStress.eof())
    {
      if (elementNoLocal != data_.functionSpace()->nElementsLocal()-1)
      {
        LOG(ERROR) << "Only load data up to element no " << elementNoLocal;
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
  this->data_.functionSpace()->geometryField().setValuesWithoutGhosts(geometryValues);
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
//! the transfer is done by the solution_vector_mapping class

typename QuasiStaticNonlinearElasticitySolverFebio::OutputConnectorDataType
QuasiStaticNonlinearElasticitySolverFebio::
getOutputConnectorData()
{
  // connect activation (input)
  ElasticitySolverOutputConnectorDataType<FieldVariableType> outputConnectorData;
  outputConnectorData.activation = this->data_.activation();

  // connect geometry (output) already done in loadFebioOutputFile


  return outputConnectorData;
}

//! output the given data for debugging

std::string QuasiStaticNonlinearElasticitySolverFebio::
getString(typename QuasiStaticNonlinearElasticitySolverFebio::OutputConnectorDataType &data)
{
  std::stringstream s;
  //s << "<QuasiStaticNonlinearElasticitySolverFebio:" << *data.activation() << ">";
  return s.str();
}

} // namespace TimeSteppingScheme
