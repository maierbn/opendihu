#include <control/computation.h>

#include <iostream>
#include "easylogging++.h"

Computation::Computation(DihuContext &dihuContext, Runnable &runnable)
{
  LOG(DEBUG)<<"start computation"<<std::endl;
  runnable.run();
}

Computation::~Computation()
{
}

void Computation::run()
{

}
