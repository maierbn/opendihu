#pragma once

#include <Python.h>  // has to be the first included header

#include <vector>

#include "control/dihu_context.h"
#include "output_writer/manager.h"
#include "basis_on_mesh/basis_on_mesh.h"

/** todo 
 */
template <int nStates>
class CellmlAdapterBase
{
public:

  ///! constructor
  CellmlAdapterBase(DihuContext context);

  ///! destructor
  ~CellmlAdapterBase();

  ///! return the compile-time constant number of state variables of one instance that will be integrated
  static constexpr int nComponents();

  ///! load model
  void initialize();
  
  ///! set initial values as given in python config
  bool setInitialValues(Vec& initialValues);

  //! return false because the object is independent of mesh type
  bool knowsMeshType();

  //! return the mesh
  std::shared_ptr<Mesh::Mesh> mesh();

  //! get number of instances, number of intermediates and number of parameters
  void getNumbers(int &nInstances, int &nIntermediates, int &nParameters);

  typedef BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<1>, BasisFunction::LagrangeOfOrder<>> BasisOnMesh;   ///< BasisOnMesh type

protected:

  //! scan the given cellml source file for initial values that are given by dummy assignments
  bool scanInitialValues(std::string sourceFilename, std::vector<double> &statesInitialValues);

  const DihuContext context_;    ///< object that contains the python config for the current context and the global singletons meshManager and solverManager
  PyObject *specificSettings_;    ///< python object containing the value of the python config dict with corresponding key
  OutputWriter::Manager outputWriterManager_; ///< manager object holding all output writer

  std::shared_ptr<Mesh::Mesh> mesh_;    ///< a mesh, there are as many instances of the same CellML problem as there are nodes in the mesh

  //int nStates_;           ///< number of states in one instance of the CellML problem (template parameter)
  int nInstances_;        ///< number of instances of the CellML problem, equals number of mesh nodes
  int nParameters_;       ///< number of parameters (=CellML name "known") in one instance of the CellML problem
  int nIntermediates_;    ///< number of intermediate values (=CellML name "wanted") in one instance of the CellML problem
   
  //std::vector<double> states_;    ///< vector of states, that are computed by rhsRoutine
  //std::vector<double> rates_;     ///< vector of rates, that are computed by rhsRoutine
  std::vector<double> parameters_; ///< vector of values that will be provided to CellML by the code, given by python config, CellML name: known
  std::vector<double> intermediates_;    ///< vector of intermediate values in DAE system. These can be computed directly from the actual states at any time. Gets computed by rhsRoutine from states, together with rates. OpenCMISS name is intermediate, CellML name: wanted

  std::string sourceFilename_; ///<file name of provided CellML right hand side routine
};

#include "cellml/00_cellml_adapter_base.tpp"