#pragma once

#include "time_stepping_scheme/time_stepping_scheme.h"
#include "output_writer/manager.h"

namespace OperatorSplitting
{

class OperatorSplitting : 
  public ::TimeSteppingScheme::TimeSteppingScheme
{
public:
  //! constructor
  OperatorSplitting(DihuContext context);
 
  virtual ~OperatorSplitting() {}

protected:
 
  OutputWriter::Manager outputWriterManager_; ///< manager object holding all output writer
  
};

}  // namespace
