#pragma once

#include <Python.h>  // has to be the first included header

#include "data_management/output_surface/output_surface.h"
#include "data_management/output_surface/convert_output_field_variables.h"

namespace Data
{

/**  The datastructures used for OutputSurface objects, these hold 2D versions of the 3D field variables
  */
template<typename Data3D>
class OutputSurface :
  public Data<typename ConvertFieldVariablesForOutputWriter<typename Data3D::FieldVariablesForOutputWriter>::FunctionSpaceFirstFieldVariable> // function space of data object is the function space of the first converted 2D field variable
{
public:

  using FirstFieldVariable = typename ConvertFieldVariablesForOutputWriter<typename Data3D::FieldVariablesForOutputWriter>::FirstFieldVariable;
  using SecondFieldVariable = typename ConvertFieldVariablesForOutputWriter<typename Data3D::FieldVariablesForOutputWriter>::SecondFieldVariable;
  using FunctionSpaceFirstFieldVariable = typename ConvertFieldVariablesForOutputWriter<typename Data3D::FieldVariablesForOutputWriter>::FunctionSpaceFirstFieldVariable;

  //! constructor
  OutputSurface(DihuContext context);

  //! initialize
  void initialize();

  //! store a pointer to the 3d data object
  void setData(Data3D &data3d);

  //! print all stored data to stdout
  void print();

  //! get if the own rank hold part of the 2D surface field variable data and thus should call the output writer for output
  bool ownRankInvolvedInOutput();

  //! get all function spaces of the faces that are extracted
  void getFunctionSpaces(std::vector<std::shared_ptr<FunctionSpaceFirstFieldVariable>> &functionSpaces);

  //! field variables that will be output by outputWriters
  typedef typename ConvertFieldVariablesForOutputWriter<typename Data3D::FieldVariablesForOutputWriter>::type FieldVariablesForOutputWriter;

  //! get pointers to all field variables that can be written by output writers
  FieldVariablesForOutputWriter getFieldVariablesForOutputWriter();

private:

  std::vector<Mesh::face_t> faces_;     ///< one of Mesh::face_t::face2Minus and Mesh::face_t::face2Plus, for which face to extract surface (other faces are not supported)
  std::shared_ptr<Data3D> data3d_;   ///< other data object that contains the 3D field variables
  bool ownRankInvolvedInOutput_;     ///< if the own rank hold part of the 2D surface field variable data and thus should call the output writer for output

  FieldVariablesForOutputWriter outputFieldVariables2D_;   ///< the surface field variables

  //! initializes the vectors with size
  void createPetscObjects() override;

};

} // namespace Data

#include "data_management/output_surface/output_surface.tpp"
