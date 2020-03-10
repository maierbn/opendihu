#include "data_management/output_surface/output_surface.h"

namespace Data
{

template<typename Data3D>
OutputSurface<Data3D>::
OutputSurface(DihuContext context) :
  Data<typename ConvertFieldVariablesForOutputWriter<typename Data3D::FieldVariablesForOutputWriter>::FunctionSpaceFirstFieldVariable>(context),
  ownRankInvolvedInOutput_(true)
{
  // parse the face from which the surface will be taken
  //std::string faceStr = this->context_.getPythonConfig().getOptionString("face","2-");

  PythonConfig specificSettings = this->context_.getPythonConfig();

  std::vector<std::string> faceStrings;
  specificSettings.template getOptionVector<std::string>("face", faceStrings);

  for (std::string faceString : faceStrings)
  {
    Mesh::face_t face = Mesh::parseFace(faceString);
    faces_.push_back(face);
    LOG(DEBUG) << "OutputSurface: parsed face \"" << Mesh::getString(face) << "\".";
  }
}

template<typename Data3D>
void OutputSurface<Data3D>::
initialize()
{
  // convert initially, this creates all 2D field variables
  typename Data3D::FieldVariablesForOutputWriter outputFieldVariables3D = this->data3d_->getFieldVariablesForOutputWriter();

  ConvertFieldVariablesForOutputWriter<typename Data3D::FieldVariablesForOutputWriter>::convert(outputFieldVariables3D, outputFieldVariables2D_, faces_, ownRankInvolvedInOutput_);

  // set function space of data3d as the function space of the first 2D field variable
  this->functionSpace_ = ConvertFieldVariablesForOutputWriter<typename Data3D::FieldVariablesForOutputWriter>::
    getFunctionSpaceFirstFieldVariable(outputFieldVariables2D_);

}

template<typename Data3D>
void OutputSurface<Data3D>::
setData(Data3D &data3d)
{
  this->data3d_ = std::make_shared<Data3D>(data3d);
}

template<typename Data3D>
void OutputSurface<Data3D>::
createPetscObjects()
{
}

template<typename Data3D>
void OutputSurface<Data3D>::
print()
{
}

template<typename Data3D>
bool OutputSurface<Data3D>::ownRankInvolvedInOutput()
{
  return ownRankInvolvedInOutput_;
}

template<typename Data3D>
void OutputSurface<Data3D>::
getFunctionSpaces(std::vector<std::shared_ptr<typename OutputSurface<Data3D>::FunctionSpaceFirstFieldVariable>> &functionSpaces)
{
  ConvertFieldVariablesForOutputWriter<typename Data3D::FieldVariablesForOutputWriter>::
    getFunctionSpacesFirstFieldVariable(outputFieldVariables2D_, functionSpaces);
}

template<typename Data3D>
typename OutputSurface<Data3D>::FieldVariablesForOutputWriter OutputSurface<Data3D>::
getFieldVariablesForOutputWriter()
{
  typename Data3D::FieldVariablesForOutputWriter outputFieldVariables3D = this->data3d_->getFieldVariablesForOutputWriter();

  // create surface field variables for 3D field variables (the template magic happens in convert_output_field_variables.h)
  ConvertFieldVariablesForOutputWriter<typename Data3D::FieldVariablesForOutputWriter>::convert(outputFieldVariables3D, outputFieldVariables2D_, faces_, ownRankInvolvedInOutput_);

  //LOG(INFO) << StringUtility::demangle(typeid(typename Data3D::FieldVariablesForOutputWriter).name());

  return outputFieldVariables2D_;
}

} // namespace
