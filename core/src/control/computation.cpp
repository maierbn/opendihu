#include <control/computation.h>

#include <iostream>
#include "easylogging++.h"

Computation::Computation(DihuContext &dihuContext, Runnable &runnable)
{
  LOG(TRACE)<<"Computation constructor: start computation by runnable.run()"<<std::endl;
  runnable.run();
}

Computation::~Computation()
{
}

void Computation::run()
{

}
