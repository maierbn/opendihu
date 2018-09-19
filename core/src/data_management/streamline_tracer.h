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
template<typename FunctionSpaceType, typename BaseDataType>
class StreamlineTracer :
  public Data<FunctionSpaceType>
{
public:

  //! constructor
  StreamlineTracer(DihuContext context);

  //! destructur
  ~StreamlineTracer();

  //! return a reference to the solution vector, the PETSc Vec can be obtained via fieldVariable->valuesGlobal()
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> gradient();

  //! print all stored data to stdout
  virtual void print();

  //! set the data object that holds the field frorm which stream lines are generated
  void setBaseData(std::shared_ptr<BaseDataType> baseData);

  //! return the number of degrees of freedom per mesh node
  static constexpr int getNDofsPerNode();

  //! return the total number of degrees of freedom, this can be a multiple of the number of nodes of the mesh
  virtual dof_no_t nNodesLocalWithGhosts();
  
  //! return the total number of degrees of freedom, this can be a multiple of the number of nodes of the mesh
  virtual dof_no_t nNodesLocalWithoutGhosts();

  //! create a fibre mesh from the given node positions, store it in mesh manager and store a pointer to the geometry field in fibreGeometry
  void createFibreMesh(const std::vector<Vec3> &nodePositions);
  
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>,BasisFunction::LagrangeOfOrder<1>> FunctionSpaceFibre;
  typedef FieldVariable::FieldVariable<FunctionSpaceFibre,3> FieldVariableFibreGeometry;
  
  //! field variables that will be output by outputWriters
  typename BaseDataType::OutputFieldVariables dummy;
  typedef decltype(std::tuple_cat(dummy, std::tuple<
    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>>,  // gradient field
    std::vector<std::shared_ptr<FieldVariableFibreGeometry>>   // geometry fields of meshes
  >())) OutputFieldVariables;

  //! get pointers to all field variables that can be written by output writers
  OutputFieldVariables getOutputFieldVariables();

protected:

  //! initializes the vectors with size
  virtual void createPetscObjects();

  std::shared_ptr<BaseDataType> baseData_;    ///< the data object that holds the field frorm which stream lines are generated
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> gradient_;    ///< the gradient field of the solution field variable

  std::vector<std::shared_ptr<FieldVariableFibreGeometry>> fibreGeometry_;   ///< geometry fields of fibres
  
  int fibreNo_ = 0; ///< counter for names of generated fibres
};

} // namespace Data

#include "data_management/streamline_tracer.tpp"
