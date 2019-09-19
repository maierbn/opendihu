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
}


void QuasiStaticNonlinearElasticitySolverFebio::
advanceTimeSpan()
{
  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  createFebioInputFile();



  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);

  // write current output values
  this->outputWriterManager_.writeOutput(this->data_, 0, endTime_);
}


void QuasiStaticNonlinearElasticitySolverFebio::
createFebioInputFile()
{
  LOG(DEBUG) << "createFebioInputFile";

  std::ofstream file("febio_input.feb");

  file << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>" << std::endl
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

    activationValue *= 1e-5;

    file << "\t\t<material id=\"" << elementNoLocal+1 << "\" name=\"Material for elmement " << elementNoLocal << "\" type=\"muscle material\">" << std::endl
      << "\t\t\t<g1>500</g1>" << std::endl
      << "\t\t\t<g2>500</g2>" << std::endl
      << "\t\t\t<p1>0.05</p1>" << std::endl
      << "\t\t\t<p2>6.6</p2>" << std::endl
      << "\t\t\t<smax>3e5</smax>" << std::endl
      << "\t\t\t<Lofl>1.07</Lofl>" << std::endl
      << "\t\t\t<lam_max>1.4</lam_max>" << std::endl
      << "\t\t\t<k>1e6</k>" << std::endl
      << "\t\t\t<fiber type=\"local\">1, 5</fiber>    <!-- fiber direction is from local node 1 to 5 in every element -->" << std::endl
      << "\t\t\t<activation lc=\"1\">" << activationValue << "</activation>" << std::endl
      << "\t\t</material>" << std::endl;
  }

  file << "\t</Material>" << std::endl
    << "\t<Geometry>" << std::endl
    << "\t\t<Nodes>" << std::endl;

  // loop over nodes
  for (node_no_t nodeNoLocal = 0; nodeNoLocal < data_.functionSpace()->nNodesLocalWithoutGhosts(); nodeNoLocal++)
  {
    dof_no_t dofNoLocal = nodeNoLocal;
    Vec3 position = data_.functionSpace()->geometryField().getValue(dofNoLocal);
    file << "\t\t\t<node id=\"" << nodeNoLocal+1 << "\">" << position[0] << "," << position[1] << "," << position[2]
      << "</node>" << std::endl;
  }

  file << "\t\t</Nodes>" << std::endl
    << "\t\t<Elements>" << std::endl;

  // loop over elements
  for (element_no_t elementNoLocal = 0; elementNoLocal < data_.functionSpace()->nElementsLocal(); elementNoLocal++)
  {
    std::array<dof_no_t,FunctionSpace::nNodesPerElement()> elementNodeNos = data_.functionSpace()->getElementNodeNos(elementNoLocal);

    file << "\t\t\t<hex8 id=\"" << elementNoLocal+1 << "\" mat=\"" << elementNoLocal+1 << "\"> ";
    file << elementNodeNos[0]+1 << ", " << elementNodeNos[2]+1 << ", " << elementNodeNos[8]+1 << ", " << elementNodeNos[6]+1 << ", "
      << elementNodeNos[18]+1 << ", " << elementNodeNos[20]+1 << ", " << elementNodeNos[26]+1 << ", " << elementNodeNos[24]+1 << ", "
      << elementNodeNos[1]+1 << ", " << elementNodeNos[5]+1 << ", " << elementNodeNos[7]+1 << ", " << elementNodeNos[3]+1 << ", "
      << elementNodeNos[19]+1 << ", " << elementNodeNos[23]+1 << ", " << elementNodeNos[25]+1 << ", " << elementNodeNos[21]+1 << ", "
      << elementNodeNos[9]+1 << ", " << elementNodeNos[11]+1 << ", " << elementNodeNos[17]+1 << ", " << elementNodeNos[15]+1;
      /*<< elementNodeNos[10] << ", " << elementNodeNos[14] << ", " << elementNodeNos[16] << ", " << elementNodeNos[12] << ", "
      << elementNodeNos[18] << ", " << elementNodeNos[13];*/
    /*for (int i = 0; i < FunctionSpace::nNodesPerElement(); i++)
    {
      if (i != 0)
        file << ", ";
      file << elementNodeNos[i];
    }*/
    file << "</hex8>" << std::endl;
  }

  file << "\t\t</Elements>" << std::endl
    << "\t</Geometry>" << std::endl
    << "\t<Boundary>" << std::endl
    << "\t\t<fix>" << std::endl;

  // fix bottom dofs
  int nNodesX = FunctionSpace::FunctionSpaceBaseDim<1,typename FunctionSpace::BasisFunction>::averageNNodesPerElement() + 1;
  int nNodesY = nNodesX;

  // fix x direction for left row
  for (int j = 0; j < nNodesY; j++)
  {
    dof_no_t dofNoLocal = j*nNodesX + nNodesX-1;
    file << "\t\t\t<node id=\"" << dofNoLocal+1 << "\" bc=\"x\"></node>" << std::endl;
  }

  file << "\t\t</fix>" << std::endl
    << "\t\t<fix> " << std::endl;

  // fix y direction for front row
  for (int i = 0; i < nNodesX; i++)
  {
    dof_no_t dofNoLocal = i;
    file << "\t\t\t<node id=\"" << dofNoLocal+1 << "\" bc=\"y\"></node>" << std::endl;
  }


  file << "\t\t</fix>" << std::endl
    << "\t\t<fix> " << std::endl;

  // fix z direction for all
  for (int j = 0; j < nNodesY; j++)
  {
    for (int i = 0; i < nNodesX; i++)
    {
      dof_no_t dofNoLocal = i;
      file << "\t\t\t<node id=\"" << dofNoLocal+1 << "\" bc=\"z\"></node>" << std::endl;
    }
  }

  file << "\t\t</fix>" << std::endl
    << "\t</Boundary>" << std::endl
    << "\t<LoadData>" << std::endl
    << "\t\t<loadcurve id=\"1\"> " << std::endl
    << "\t\t\t<loadpoint>0,0</loadpoint>" << std::endl
    << "\t\t\t<loadpoint>1,1e-5</loadpoint>" << std::endl
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
    << "\t\t\t<node_data file=\"febio_geometry_output.txt\" format=\"%i,%g,%g,%g\" data=\"x;y;z\"/>" << std::endl
    << "\t\t</logfile>" << std::endl
    << "\t</Output>" << std::endl
    << "</febio_spec>" << std::endl;

  file.close();
}

void QuasiStaticNonlinearElasticitySolverFebio::
loadFebioOutputFile()
{

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
  ElasticitySolverOutputConnectorDataType<FieldVariableType> outputConnectorData;
  outputConnectorData.activation = this->data_.activation();

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
