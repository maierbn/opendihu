#include "specialized_solver/solid_mechanics/quasi_static/febio/nonlinear_elasticity_solver_febio.h"

#include <Python.h>  // has to be the first included header

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "data_management/specialized_solver/multidomain.h"
#include "control/diagnostic_tool/performance_measurement.h"
#include "control/diagnostic_tool/solver_structure_visualizer.h"

namespace TimeSteppingScheme
{
#if 1
NonlinearElasticitySolverFebio::
NonlinearElasticitySolverFebio(DihuContext context, std::string solverName) :
  context_(context[solverName]), data_(context_), solverName_(solverName), initialized_(false)
{

  // get python config
  this->specificSettings_ = this->context_.getPythonConfig();

  // read in the durationLogKey
  if (specificSettings_.hasKey("durationLogKey"))
  {
    this->durationLogKey_ = specificSettings_.getOptionString("durationLogKey", "");
  }

  // load settings
  this->activationFactor_ = 0;

  // parse material parameters
  if (solverName_ == "NonlinearElasticitySolverFebio")    // if the base class solver is referenced in settings
  {
    specificSettings_.getOptionVector("materialParameters", materialParameters_);

    if (materialParameters_.size() < 3)
    {
      LOG(FATAL) << specificSettings_ << "[\"materialParameters\"]: 3 material parameters, [c0, c1, k], are required, but only " << materialParameters_.size() << " are given.";
    }
  }

  // load traction elements
  specificSettings_.getOptionVector("tractionElementNos", tractionElementNos_);
  tractionVector_ = specificSettings_.getOptionArray<double,3>("tractionVector", 0);

  // initialize output writers
  this->outputWriterManager_.initialize(this->context_, this->specificSettings_);

  problemDescription_ = "Mesh to solve";
  problemTitle_ = "Isotropic Mooney-Rivlin";

  LOG(DEBUG) << "initialized NonlinearElasticitySolverFebio";
}


void NonlinearElasticitySolverFebio::
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
  this->outputWriterManager_.writeOutput(this->data_, 1, endTime_);
}

bool NonlinearElasticitySolverFebio::
isFebioAvailable()
{
  // to see if febio2 is installed, run it with a non-existing command line argument "a"
  // this produces the output "FATAL ERROR: Invalid command line option" and exits with status 0
  int returnValue = system("febio3 a > /dev/null");

  return returnValue == EXIT_SUCCESS;
}

void NonlinearElasticitySolverFebio::
runFebio()
{
  // only run febio on rank 0
  int ownRankNo = DihuContext::ownRankNoCommWorld();

  if (ownRankNo == 0)
  {

    // run febio simulation
    int returnValue = system("rm -f febio3_input.log febio3_stress_output.txt febio3_geometry_output.txt; febio3 -i febio3_input.feb > /dev/null");

    if (returnValue != EXIT_SUCCESS)
    {
      LOG(ERROR) << "Running febio failed with error code " << returnValue;
    }

#ifndef NDEBUG
    // copy contents of log file from the febio execution to the opendihu.log file
    std::ifstream logFile("febio3_input.log");

    if (logFile.is_open())
    {
      std::string logContents((std::istreambuf_iterator<char>(logFile)), std::istreambuf_iterator<char>());
      LOG(DEBUG) << "Content of febio3 log: \n" << logContents;
    }
#endif
  }

  // barrier to wait until febio simulation is finished
  MPIUtility::handleReturnValue(MPI_Barrier(this->data_.functionSpace()->meshPartition()->mpiCommunicator()), "MPI_Barrier");
}

void NonlinearElasticitySolverFebio::
createFebioInputFile()
{
  LOG(DEBUG) << "createFebioInputFile base";

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

    fileContents << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
      << "<!--\n"
      << "Problem Description:\n"
      << "\t" << problemDescription_ << " at time " << endTime_ << ", created by opendihu\n"
      << "-->\n"
      << "<febio_spec version=\"2.5\">\n"
      << "\t<Module type=\"solid\"/>\n"
      << "\t<Control>\n"
      << "\t\t<title>" << problemTitle_ << "</title>\n"
  //    << "\t\t<restart file=\"dump.out\">1</restart>\n"
      << "\t\t<time_steps>10</time_steps>\n"
      << "\t\t<step_size>0.1</step_size>\n"
      << "\t\t<max_refs>100</max_refs> <!-- (15) Max number of stiffness reformations -->\n"
      << "\t\t<max_ups>10</max_ups>   <!-- (10) Max number of BFGS/Broyden stiffness updates --> \n"
      << "\t\t<dtol>1e-10</dtol>      <!-- (0.001) Convergence tolerance on displacement -->\n"
      << "\t\t<etol>1e-10</etol>       <!-- (0.01) Convergence tolerance on energy -->\n"
      << "\t\t<rtol>0</rtol>          <!-- Convergence tolerance on residual, 0=disabled-->\n"
      << "\t\t<lstol>1e-10</lstol>      <!-- (0.9) Convergence tolerance on line search -->\n"
      << "\t\t<analysis type=\"static\"></analysis>\n"
      << "\t\t<time_stepper>\n"
      << "\t\t\t<dtmin>0.01</dtmin>\n"
      << "\t\t\t<dtmax>0.1</dtmax>\n"
      << "\t\t\t<max_retries>5</max_retries>\n"
      << "\t\t\t<opt_iter>10</opt_iter>\n"
      << "\t\t</time_stepper>\n"
      << "\t</Control>\n"
      << "\t<Material>\n"
      << "\t\t<material id=\"1\" name=\"Material\" type=\"Mooney-Rivlin\">\n"
      << "\t\t\t<c1>" << materialParameters_[0] << "</c1>\n"
      << "\t\t\t<c2>" << materialParameters_[1] << "</c2>\n"
      << "\t\t\t<k>" << materialParameters_[2] << "</k>\n"
      << "\t\t</material>\n"
      << "\t</Material>\n"
      << "\t<Geometry>\n"
      << "\t\t<Nodes name=\"febio3_input\">\n";
    // loop over nodes
    for (global_no_t nodeNoGlobalPetsc = 0; nodeNoGlobalPetsc < data_.functionSpace()->nNodesGlobal(); nodeNoGlobalPetsc++)
    {
      double nodePositionX = nodePositionValuesGlobal[3*nodeNoGlobalPetsc + 0];
      double nodePositionY = nodePositionValuesGlobal[3*nodeNoGlobalPetsc + 1];
      double nodePositionZ = nodePositionValuesGlobal[3*nodeNoGlobalPetsc + 2];

      fileContents << "\t\t\t<node id=\"" << nodeNoGlobalPetsc+1 << "\">"
        << nodePositionX << "," << nodePositionY << "," << nodePositionZ << "</node>\n";
    }

    fileContents << "\t\t</Nodes>\n"
      << "\t\t<Elements type=\"hex8\" mat=\"1\" name=\"Part1\">\n";

    // loop over global elements and insert node nos of elements
    for (global_no_t elementNoGlobalPetsc = 0; elementNoGlobalPetsc < data_.functionSpace()->nElementsGlobal(); elementNoGlobalPetsc++)
    {
      nodeNosGlobal[elementNoGlobalPetsc*8 + 0];

      // elements types: https://help.febio.org/FEBio/FEBio_um_2_8/FEBio_um_2-8-3.8.2.1.html
      fileContents << "\t\t\t<elem id=\"" << elementNoGlobalPetsc+1 << "\"> ";
      fileContents << nodeNosGlobal[elementNoGlobalPetsc*8 + 0] << ", " << nodeNosGlobal[elementNoGlobalPetsc*8 + 1]
        << ", " << nodeNosGlobal[elementNoGlobalPetsc*8 + 2] << ", " << nodeNosGlobal[elementNoGlobalPetsc*8 + 3]
        << ", " << nodeNosGlobal[elementNoGlobalPetsc*8 + 4] << ", " << nodeNosGlobal[elementNoGlobalPetsc*8 + 5]
        << ", " << nodeNosGlobal[elementNoGlobalPetsc*8 + 6] << ", " << nodeNosGlobal[elementNoGlobalPetsc*8 + 7];
      fileContents << "</elem>\n";
    }

    fileContents << "\t\t</Elements>\n"
      << "\t\t<NodeSet name=\"FixedDisplacementX\">\n";

    // fix bottom dofs
    int nNodesX = this->data_.functionSpace()->nNodesGlobal(0);
    int nNodesY = this->data_.functionSpace()->nNodesGlobal(1);
    //int nNodesZ = this->data_.functionSpace()->nNodesGlobal(2);
    //int nElementsXY = (nNodesX-1)*(nNodesY-1);

    enum {fixAll, fixFloating} dirichletBoundaryConditionsMode = fixFloating;
    std::string dirichletBoundaryConditionsModeString = this->specificSettings_.getOptionString("dirichletBoundaryConditionsMode", "fix_all");

    if (dirichletBoundaryConditionsModeString == "fix_all")
      dirichletBoundaryConditionsMode = fixAll;

    if (dirichletBoundaryConditionsMode == fixAll)
    {
      // fix x direction for all
      for (int j = 0; j < nNodesY; j++)
      {
        for (int i = 0; i < nNodesX; i++)
        {
          dof_no_t dofNoLocal = j*nNodesX + i;
          fileContents << "\t\t\t<node id=\"" << dofNoLocal+1 << "\" />\n";
        }
      }
    }
    else
    {
      // fix x direction for left row
      for (int j = 0; j < nNodesY; j++)
      {
        dof_no_t dofNoLocal = j*nNodesX;
        fileContents << "\t\t\t<node id=\"" << dofNoLocal+1 << "\" />\n";
      }
    }

    fileContents << "\t\t</NodeSet>\n"
      << "\t\t<NodeSet name=\"FixedDisplacementY\">\n";

    if (dirichletBoundaryConditionsMode == fixAll)
    {
      // fix x direction for all
      for (int j = 0; j < nNodesY; j++)
      {
        for (int i = 0; i < nNodesX; i++)
        {
          dof_no_t dofNoLocal = j*nNodesX + i;
          fileContents << "\t\t\t<node id=\"" << dofNoLocal+1 << "\" />\n";
        }
      }
    }
    else
    {
      // fix y direction for front row
      for (int i = 0; i < nNodesX; i++)
      {
        dof_no_t dofNoLocal = i;
        fileContents << "\t\t\t<node id=\"" << dofNoLocal+1 << "\" />\n";
      }
    }

    fileContents << "\t\t</NodeSet>\n"
      << "\t\t<NodeSet name=\"FixedDisplacementZ\">\n";

    if (dirichletBoundaryConditionsMode == fixAll)
    {
      // fix x direction for all
      for (int j = 0; j < nNodesY; j++)
      {
        for (int i = 0; i < nNodesX; i++)
        {
          dof_no_t dofNoLocal = j*nNodesX + i;
          fileContents << "\t\t\t<node id=\"" << dofNoLocal+1 << "\" />\n";
        }
      }
    }
    else
    {
      // fix z direction for all
      for (int j = 0; j < nNodesY; j++)
      {
        for (int i = 0; i < nNodesX; i++)
        {
          dof_no_t dofNoLocal = j*nNodesX + i;
          fileContents << "\t\t\t<node id=\"" << dofNoLocal+1 << "\" />\n";
        }
      }
    }

    fileContents << "\t\t</NodeSet>\n"
      << "\t\t<Surface name=\"SurfaceTraction1\">\n";

    // iterate over surface traction elements
    for (element_no_t elementNoLocal : tractionElementNos_)
    {
      fileContents << "\t\t\t<quad4 id=\"1\">";

      // loop over 4 top nodes of element
      // febio  opendihu
      //  3 2     6 7
      //  0 1     4 5
      std::array<int,4> dofIndices = {4,5,7,6};


      for (int dofIndex : dofIndices)
      {
        // get global dof nos of the current element
        dof_no_t dofNoLocal = this->data_.functionSpace()->getDofNo(elementNoLocal, dofIndex);
        global_no_t dofNoGlobalPetsc = this->data_.functionSpace()->meshPartition()->getDofNoGlobalPetsc(dofNoLocal);

        if (dofIndex > 4)
          fileContents << ", ";
        fileContents << (dofNoGlobalPetsc+1);
      }
      fileContents << "</quad4>\n";
    }


    fileContents << "\t\t</Surface>\n"
      << "\t</Geometry>\n"
      << "\t<Boundary>\n"
      << "\t\t<fix bc=\"x\" node_set=\"FixedDisplacementX\"/>\n"
      << "\t\t<fix bc=\"y\" node_set=\"FixedDisplacementY\"/>\n"
      << "\t\t<fix bc=\"z\" node_set=\"FixedDisplacementZ\"/>\n"
      << "\t</Boundary>\n"
      << "\t<Loads>\n";
    fileContents
      << "\t\t<surface_load type=\"traction\" surface=\"SurfaceTraction1\">\n"
      << "\t\t\t<scale lc=\"2\">1</scale>\n"
      << "\t\t\t<traction>" << tractionVector_[0] << "," << tractionVector_[1] << "," << tractionVector_[2] << "</traction>\n"
      << "\t\t</surface_load>\n"
      << "\t</Loads>\n";

    fileContents << "\t<LoadData>\n"
      << "\t\t<loadcurve id=\"1\"> <!-- for activation -->\n"
      << "\t\t\t<loadpoint>0,0</loadpoint>\n"
      << "\t\t\t<loadpoint>1," << activationFactor_ << "</loadpoint>\n"
      << "\t\t</loadcurve>\n"
      << "\t\t<loadcurve id=\"2\"> <!-- for external load (surface traction) that stretches the muscle -->\n"
      << "\t\t\t<loadpoint>0,0</loadpoint>\n"
      << "\t\t\t<loadpoint>1,1</loadpoint>\n"
      << "\t\t</loadcurve>\n"
      << "\t</LoadData>\n"
      << "\t<Output>\n"
      << "\t\t<plotfile type=\"febio\">\n"
      << "\t\t\t<var type=\"fiber vector\"/>\n"
      << "\t\t\t<var type=\"displacement\"/>\n"
      << "\t\t\t<var type=\"fiber stretch\"/>\n"
      << "\t\t\t<var type=\"Lagrange strain\"/>\n"
      << "\t\t\t<var type=\"parameter\"/>\n"
      << "\t\t\t<var type=\"relative volume\"/>\n"
      << "\t\t\t<var type=\"stress\"/>\n"
      << "\t\t\t<var type=\"strain energy density\"/>\n"
      << "\t\t</plotfile>\n"
      << "\t\t<logfile>\n"

      // available variables: https://help.febio.org/FEBio/FEBio_um_2_9/index.html Sec. 3.17.1.2 and 3.17.1.3
      << "\t\t\t<node_data file=\"febio3_geometry_output.txt\" format=\"%i,%g,%g,%g,%g,%g,%g,%g,%g,%g\" data=\"x;y;z;ux;uy;uz;Rx;Ry;Rz\"/>\n"
      << "\t\t\t<element_data file=\"febio3_stress_output.txt\" format=\"%i,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\" data=\"sx;sy;sz;sxy;syz;sxz;Ex;Ey;Ez;Exy;Eyz;Exz;J;Fxx;Fxy;Fxz;Fyx;Fyy;Fyz;Fzx;Fzy;Fzz\"/>\n"
      << "\t\t</logfile>\n"
      << "\t</Output>\n"
      << "</febio_spec>\n";

    std::ofstream file("febio3_input.feb");

    if (!file.is_open())
    {
      LOG(ERROR) << "Could not write to file \"febio3_input.feb\".";
    }

    file << fileContents.str();
    file.close();
  }
}

void NonlinearElasticitySolverFebio::
loadFebioOutputFile()
{
  LOG(TRACE) << "loadFebioOutputFile";

  std::ifstream fileGeometry;
  fileGeometry.open("febio3_geometry_output.txt", std::ios::in | std::ios::binary);

  if (!fileGeometry.is_open())
  {
    LOG(WARNING) << "Could not read output file \"febio3_geometry_output.txt\" that should have been created by FEBio.";
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
  fileStress.open("febio3_stress_output.txt", std::ios::in | std::ios::binary);

  if (!fileStress.is_open())
  {
    LOG(WARNING) << "Could not read output file \"febio3_stress_output.txt\" that should have been created by FEBio.";
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

  this->data_.pk2Stress()->zeroEntries();
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

      // sx;sy;sz;sxy;syz;sxz;Ex;Ey;Ez;Exy;Eyz;Exz;J;Fxx,Fxy,Fxz;Fyx;Fyy;Fyz;Fzx;Fzy;Fzz
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
      double J   = atof(StringUtility::extractUntil(line, ",").c_str());
      double Fxx = atof(StringUtility::extractUntil(line, ",").c_str());
      double Fxy = atof(StringUtility::extractUntil(line, ",").c_str());
      double Fxz = atof(StringUtility::extractUntil(line, ",").c_str());
      double Fyx = atof(StringUtility::extractUntil(line, ",").c_str());
      double Fyy = atof(StringUtility::extractUntil(line, ",").c_str());
      double Fyz = atof(StringUtility::extractUntil(line, ",").c_str());
      double Fzx = atof(StringUtility::extractUntil(line, ",").c_str());
      double Fzy = atof(StringUtility::extractUntil(line, ",").c_str());
      double Fzz = atof(line.c_str());

      // compute 2nd Piola-Kirchhoff stress, S, from Cauchy stress, σ
      // S = J F^-1 σ F^-T
      Tensor2<3> cauchyStress{Vec3{sx,sxy,sxz}, Vec3{sxy, sy, syz}, Vec3{sxz, syz, sz}};
      Tensor2<3> deformationGradient{Vec3{Fxx, Fyx, Fzx}, Vec3{Fxy, Fyy, Fzy}, Vec3{Fxz, Fyz, Fzz}};
      double determinant = 0;
      Tensor2<3> inverseDeformationGradient = MathUtility::computeInverse(deformationGradient, determinant);

      Tensor2<3> deformationGradientCofactor = MathUtility::computeCofactorMatrix<double>(deformationGradient);  // cof(M) = det(M) * M^{-T}
      Tensor2<3> pk2Stress = inverseDeformationGradient * cauchyStress * deformationGradientCofactor;

      LOG(DEBUG) << "local element " << elementNoLocal << " of " << this->data_.functionSpace()->nElementsLocal() << ", " << FunctionSpace::nNodesPerElement() << " elementNodeNos: " << elementNodeNos;

      LOG(DEBUG) << "F: " << deformationGradient << ", J: " << J << "=" << determinant;
      LOG(DEBUG) << "pk2Stress = " << pk2Stress << " = " << inverseDeformationGradient << "*" << cauchyStress << "*" << deformationGradientCofactor;

      for (node_no_t elementalNodeNo = 0; elementalNodeNo < FunctionSpace::nNodesPerElement(); elementalNodeNo++)
      {
        node_no_t nodeNoLocal = elementNodeNos[elementalNodeNo];

        if (nodeNoLocal >= this->data_.functionSpace()->nNodesLocalWithoutGhosts())
          continue;

        this->data_.pk2Stress()->setValue(nodeNoLocal, std::array<double,6>{pk2Stress[0][0],pk2Stress[1][1],pk2Stress[2][2],pk2Stress[1][0],pk2Stress[2][1],pk2Stress[2][0]}, ADD_VALUES);
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
        LOG(ERROR) << "Loaded " << nElementsLoaded << " local elements from \"febio3_stress_output.txt\", but " << this->data_.functionSpace()->nElementsLocal() << " were expected.";
      }
      break;
    }
  }

  // divide values at nodes by number of summands
  // loop over nodes
  for (node_no_t nodeNoLocal = 0; nodeNoLocal < data_.functionSpace()->nNodesLocalWithoutGhosts(); nodeNoLocal++)
  {
    // PK2 stress
    std::array<double,6> pk2StressValues = this->data_.pk2Stress()->getValue(nodeNoLocal);
    pk2StressValues /= nSummands[nodeNoLocal];
    this->data_.pk2Stress()->setValue(nodeNoLocal, pk2StressValues, INSERT_VALUES);

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

  //this->data_.functionSpace()->geometryField().finishGhostManipulation();
  this->data_.functionSpace()->geometryField().setValuesWithoutGhosts(geometryValues);

  this->data_.functionSpace()->geometryField().zeroGhostBuffer();
  this->data_.functionSpace()->geometryField().setRepresentationGlobal();
  this->data_.functionSpace()->geometryField().startGhostManipulation();

  LOG(DEBUG) << "geometryField pointer: " << this->data_.functionSpace()->geometryField().partitionedPetscVec();
  LOG(DEBUG) << "referenceGeometry pointer: " << this->data_.referenceGeometry()->partitionedPetscVec();
}

void NonlinearElasticitySolverFebio::
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

  //std::cout << "at StimulationLogging::logStimulationBegin, ownRankNo: " << ownRankNo << ", nRanks: " << nRanks << "\n";
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

void NonlinearElasticitySolverFebio::
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

  //std::cout << "at StimulationLogging::logStimulationBegin, ownRankNo: " << ownRankNo << ", nRanks: " << nRanks << "\n";
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

void NonlinearElasticitySolverFebio::
run()
{
  // initialize everything
  initialize();

  this->advanceTimeSpan();
}


void NonlinearElasticitySolverFebio::
setTimeSpan(double startTime, double endTime)
{
  endTime_ = endTime;
}


void NonlinearElasticitySolverFebio::
initialize()
{
  LOG(DEBUG) << "initialize NonlinearElasticitySolverFebio";
  assert(this->specificSettings_.pyObject());

  // create function space
  std::shared_ptr<FunctionSpace> functionSpace = context_.meshManager()->functionSpace<FunctionSpace>(specificSettings_);

  // initialize the data object
  // store mesh in data
  data_.setFunctionSpace(functionSpace);
  data_.initialize();

  // write initial geometry
  this->outputWriterManager_.writeOutput(this->data_, 0, 0);

  // add this solver to the solvers diagram
  DihuContext::solverStructureVisualizer()->addSolver(solverName_);

  // set the outputConnectorData for the solverStructureVisualizer to appear in the solver diagram
  DihuContext::solverStructureVisualizer()->setOutputConnectorData(getOutputConnectorData());

  LOG(DEBUG) << "initialization done";
  this->initialized_ = true;
}


void NonlinearElasticitySolverFebio::reset()
{
  this->initialized_ = false;
}


typename NonlinearElasticitySolverFebio::Data &NonlinearElasticitySolverFebio::
data()
{
  return data_;
}

//! get the data that will be transferred in the operator splitting to the other term of the splitting
//! the transfer is done by the output_connector_data_transfer class

std::shared_ptr<typename NonlinearElasticitySolverFebio::OutputConnectorDataType>
NonlinearElasticitySolverFebio::
getOutputConnectorData()
{
  return this->data_.getOutputConnectorData();
}

//! output the given data for debugging

std::string NonlinearElasticitySolverFebio::
getString(std::shared_ptr<typename NonlinearElasticitySolverFebio::OutputConnectorDataType> data)
{
  std::stringstream s;
  //s << "<NonlinearElasticitySolverFebio:" << *data.activation() << ">";
  return s.str();
}
#endif    // disable whole solver, because it takes long to compile

} // namespace TimeSteppingScheme
