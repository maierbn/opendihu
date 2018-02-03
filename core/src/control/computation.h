#pragma once

#include "dihu_context.h"
#include "runnable.h"

class Computation : public Runnable
{
public:
  Computation(DihuContext &dihuContext, Runnable &runnable);
  virtual ~Computation();
  
  // run the computation
  void run();
private:
  
};
