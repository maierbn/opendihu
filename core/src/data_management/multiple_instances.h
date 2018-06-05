#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "control/types.h"
#include "mesh/mesh.h"
#include "data_management/data.h"
#include "field_variable/field_variable.h"

namespace Data
{

/**  The datastructures used for streamline tracer.
 *   BaseDataType is a Data class that provides the solution field variable for the streamline tracer to operate on.
 */
template<typename BasisOnMeshType, typename BaseTimesteppingType>
class MultipleInstances :
  public Data<BasisOnMeshType>
{
public:
  typedef typename BaseTimesteppingType::Data BaseDataType;

  //! constructor
  MultipleInstances(DihuContext context);

  //! destructur
  ~MultipleInstances();

  //! print all stored data to stdout
  virtual void print();

  //! set the data objects of the instances
  void setInstancesData(std::vector<BaseTimesteppingType> &instances);

  //! return the number of degrees of freedom per mesh node
  static constexpr int getNDofsPerNode();

  //! return the total number of degrees of freedom, this can be a multiple of the number of nodes of the mesh
  virtual dof_no_t nUnknowns();

  typedef BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<1>,BasisFunction::LagrangeOfOrder<1>> MeshFibre;
  typedef FieldVariable::FieldVariable<MeshFibre,3> FieldVariableFibreGeometry;
  
  //! field variables that will be output by outputWriters
  typedef std::tuple<std::vector<typename BaseDataType::OutputFieldVariables>> OutputFieldVariables;

  //! get pointers to all field variables that can be written by output writers
  OutputFieldVariables getOutputFieldVariables();

protected:

  //! initializes the vectors with size
  virtual void createPetscObjects();

  std::vector<std::shared_ptr<BaseDataType>> instancesData_;    ///< the data objects of all instances
};

} // namespace Data

#include "data_management/multiple_instances.tpp"
