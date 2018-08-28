#include "output_writer/python_callback/python_callback.h"

#include <Python.h>  // has to be the first included header
#include <iostream>

#include "easylogging++.h"
#include "output_writer/python_callback/python_callback_writer.h"

namespace OutputWriter
{

template<typename DataType>
void PythonCallback::write(DataType& data, int timeStepNo, double currentTime)
{
  // check if output should be written in this timestep and prepare filename
  if (!Generic::prepareWrite(data, timeStepNo, currentTime))
  {
    return;
  }
  typedef typename DataType::BasisOnMesh MeshType;
  std::shared_ptr<MeshType> mesh = std::static_pointer_cast<MeshType>(data.mesh());

  LOG(TRACE) << "PythonCallback::write timeStepNo=" <<timeStepNo<< ", currentTime=" << currentTime;

  // call implementation specific for BasisOnMesh type
  PythonCallbackWriter<typename DataType::BasisOnMesh,typename DataType::OutputFieldVariables>::
    callCallback(callback_, data.getOutputFieldVariables(), this->timeStepNo_, this->currentTime_, this->onlyNodalValues_);
}

};