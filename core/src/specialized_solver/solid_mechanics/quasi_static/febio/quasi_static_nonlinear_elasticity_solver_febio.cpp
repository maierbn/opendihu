#include "specialized_solver/solid_mechanics/quasi_static/febio/quasi_static_nonlinear_elasticity_solver_febio.h"

#include <Python.h>  // has to be the first included header

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "data_management/specialized_solver/multidomain.h"
#include "control/diagnostic_tool/performance_measurement.h"
#include "control/diagnostic_tool/solver_structure_visualizer.h"

namespace TimeSteppingScheme
{
#if 1
QuasiStaticNonlinearElasticitySolverFebio::
QuasiStaticNonlinearElasticitySolverFebio(DihuContext context) :
  NonlinearElasticitySolverFebio(context, "QuasiStaticNonlinearElasticitySolverFebio")
{
  // load settings
  this->activationFactor_ = this->specificSettings_.getOptionDouble("activationFactor", 1e-5, PythonUtility::NonNegative);
  this->force_ = this->specificSettings_.getOptionDouble("force", 100, PythonUtility::NonNegative);

  LOG(DEBUG) << "initialized QuasiStaticNonlinearElasticitySolverFebio";

  this->problemDescription_ = "Muscle mesh to solve";
  this->problemTitle_ = "Biceps";
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

    fileContents << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>" << "\n"
      << "<!--" << "\n"
      << "Problem Description:" << "\n"
      << "\t" << problemDescription_ << " at time " << endTime_ << ", created by opendihu" << "\n"
      << "-->" << "\n"
      << "<febio_spec version=\"1.2\">" << "\n"
      << "\t<Control>" << "\n"
      << "\t\t<title>" << problemTitle_ << "</title>" << "\n"
  //    << "\t\t<restart file=\"dump.out\">1</restart>" << "\n"
      << "\t\t<time_steps>10</time_steps>" << "\n"
      << "\t\t<step_size>0.1</step_size>" << "\n"
      << "\t\t<max_refs>15</max_refs> <!-- Max number of stiffness reformations -->" << "\n"
      << "\t\t<max_ups>10</max_ups>   <!-- Max number of BFGS/Broyden stiffness updates --> " << "\n"
      << "\t\t<dtol>0.001</dtol>      <!-- Convergence tolerance on displacement -->" << "\n"
      << "\t\t<etol>0.01</etol>       <!-- Convergence tolerance on energy -->" << "\n"
      << "\t\t<rtol>0</rtol>          <!-- Convergence tolerance on residual, 0=disabled-->" << "\n"
      << "\t\t<lstol>0.9</lstol>      <!-- Convergence tolerance on line search -->" << "\n"
      << "\t\t<analysis type=\"dynamic\"></analysis>" << "\n"
      << "\t\t<time_stepper>" << "\n"
      << "\t\t\t<dtmin>0.01</dtmin>" << "\n"
      << "\t\t\t<dtmax>0.1</dtmax>" << "\n"
      << "\t\t\t<max_retries>5</max_retries>" << "\n"
      << "\t\t\t<opt_iter>10</opt_iter>" << "\n"
      << "\t\t</time_stepper>" << "\n"
      << "\t</Control>" << "\n"
      << "\t<Material>" << "\n";

    // loop over global elements and add one material per element with the activation value of this element
    for (global_no_t elementNoGlobalPetsc = 0; elementNoGlobalPetsc < data_.functionSpace()->nElementsGlobal(); elementNoGlobalPetsc++)
    {
      double activationValue = activationValuesGlobal[elementNoGlobalPetsc];

      fileContents << "\t\t<material id=\"" << elementNoGlobalPetsc+1 << "\" name=\"Material for element " << elementNoGlobalPetsc << "\" type=\"muscle material\">" << "\n"
        << "\t\t\t<g1>500</g1>" << "\n"
        << "\t\t\t<g2>500</g2>" << "\n"
        << "\t\t\t<p1>0.05</p1>" << "\n"
        << "\t\t\t<p2>6.6</p2>" << "\n"
        << "\t\t\t<smax>3e5</smax>" << "\n"
        << "\t\t\t<Lofl>1.07</Lofl>" << "\n"
        << "\t\t\t<lam_max>1.4</lam_max>" << "\n"
        << "\t\t\t<k>1e6</k>" << "\n"
        << "\t\t\t<fiber type=\"local\">1, 5</fiber>    <!-- fiber direction is from local node 1 to 5 in every element, i.e. upwards -->" << "\n"
        << "\t\t\t<activation lc=\"1\">" << activationValue << "</activation>" << "\n"
        << "\t\t</material>" << "\n";
    }

    fileContents << "\t</Material>" << "\n"
      << "\t<Geometry>" << "\n"
      << "\t\t<Nodes>" << "\n";

    // loop over nodes
    for (global_no_t nodeNoGlobalPetsc = 0; nodeNoGlobalPetsc < data_.functionSpace()->nNodesGlobal(); nodeNoGlobalPetsc++)
    {
      double nodePositionX = nodePositionValuesGlobal[3*nodeNoGlobalPetsc + 0];
      double nodePositionY = nodePositionValuesGlobal[3*nodeNoGlobalPetsc + 1];
      double nodePositionZ = nodePositionValuesGlobal[3*nodeNoGlobalPetsc + 2];

      fileContents << "\t\t\t<node id=\"" << nodeNoGlobalPetsc+1 << "\">"
        << nodePositionX << "," << nodePositionY << "," << nodePositionZ << "</node>" << "\n";
    }

    fileContents << "\t\t</Nodes>" << "\n"
      << "\t\t<Elements>" << "\n";

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
      fileContents << "</hex8>" << "\n";
    }

    fileContents << "\t\t</Elements>" << "\n"
      << "\t</Geometry>" << "\n"
      << "\t<Boundary>" << "\n"
      << "\t\t<fix>" << "\n";

    // fix bottom dofs
    int nNodesX = this->data_.functionSpace()->nNodesGlobal(0);
    int nNodesY = this->data_.functionSpace()->nNodesGlobal(1);
    int nNodesZ = this->data_.functionSpace()->nNodesGlobal(2);

    // fix x direction for left row
    for (int j = 0; j < nNodesY; j++)
    {
      dof_no_t dofNoLocal = j*nNodesX;
      fileContents << "\t\t\t<node id=\"" << dofNoLocal+1 << "\" bc=\"x\"></node>" << "\n";
    }

    fileContents << "\t\t</fix>" << "\n"
      << "\t\t<fix> " << "\n";

    // fix y direction for front row
    for (int i = 0; i < nNodesX; i++)
    {
      dof_no_t dofNoLocal = i;
      fileContents << "\t\t\t<node id=\"" << dofNoLocal+1 << "\" bc=\"y\"></node>" << "\n";
    }


    fileContents << "\t\t</fix>" << "\n"
      << "\t\t<fix> " << "\n";

    // fix z direction for all
    for (int j = 0; j < nNodesY; j++)
    {
      for (int i = 0; i < nNodesX; i++)
      {
        dof_no_t dofNoLocal = j*nNodesX + i;
        fileContents << "\t\t\t<node id=\"" << dofNoLocal+1 << "\" bc=\"z\"></node>" << "\n";
      }
    }

    fileContents << "\t\t</fix>" << "\n"
      << "\t</Boundary>" << "\n";

    // force that stretches the muscle
    fileContents << "\t<Loads>" << "\n"
      << "\t\t<force>" << "\n";

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
          << value << "</node>" << "\n";
      }
    }


    fileContents << "\t\t</force>" << "\n"
      << "\t</Loads>" << "\n";

    fileContents << "\t<LoadData>" << "\n"
      << "\t\t<loadcurve id=\"1\"> <!-- for activation -->" << "\n"
      << "\t\t\t<loadpoint>0,0</loadpoint>" << "\n"
      << "\t\t\t<loadpoint>1," << activationFactor_ << "</loadpoint>" << "\n"
      << "\t\t</loadcurve>" << "\n"
      << "\t\t<loadcurve id=\"2\"> <!-- for external load that stretches the muscle -->" << "\n"
      << "\t\t\t<loadpoint>0,0</loadpoint>" << "\n"
      << "\t\t\t<loadpoint>1," << force_ << "</loadpoint>" << "\n"
      << "\t\t</loadcurve>" << "\n"
      << "\t</LoadData>" << "\n"
      << "\t<Output>" << "\n"
      << "\t\t<plotfile type=\"febio\">" << "\n"
      << "\t\t\t<var type=\"fiber vector\"/>" << "\n"
      << "\t\t\t<var type=\"displacement\"/>" << "\n"
      << "\t\t\t<var type=\"fiber stretch\"/>" << "\n"
      << "\t\t\t<var type=\"Lagrange strain\"/>" << "\n"
      << "\t\t\t<var type=\"parameter\"/>" << "\n"
      << "\t\t\t<var type=\"relative volume\"/>" << "\n"
      << "\t\t\t<var type=\"stress\"/>" << "\n"
      << "\t\t\t<var type=\"strain energy density\"/>" << "\n"
      << "\t\t</plotfile>" << "\n"
      << "\t\t<logfile>" << "\n"

      // available variables: https://help.febio.org/FEBio/FEBio_um_2_9/index.html Sec. 3.17.1.2 and 3.17.1.3
      << "\t\t\t<node_data file=\"febio3_geometry_output.txt\" format=\"%i,%g,%g,%g,%g,%g,%g,%g,%g,%g\" data=\"x;y;z;ux;uy;uz;Rx;Ry;Rz\"/>" << "\n"
      << "\t\t\t<element_data file=\"febio3_stress_output.txt\" format=\"%i,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\" data=\"sx;sy;sz;sxy;syz;sxz;Ex;Ey;Ez;Exy;Eyz;Exz;J;Fxx;Fxy;Fxz;Fyx;Fyy;Fyz;Fzx;Fzy;Fzz\"/>" << "\n"
      << "\t\t</logfile>" << "\n"
      << "\t</Output>" << "\n"
      << "</febio_spec>" << "\n";


    std::ofstream file("febio3_input.feb");

    if (!file.is_open())
    {
      LOG(ERROR) << "Could not write to file \"febio3_input.feb\".";
    }

    file << fileContents.str();
    file.close();
  }
}

#endif    // disable whole solver, because it takes long to compile

} // namespace TimeSteppingScheme
