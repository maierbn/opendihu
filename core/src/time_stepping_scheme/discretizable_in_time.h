#pragma once

#include <petscsys.h>

class DiscretizableInTime
{
public:
  virtual void evaluateTimesteppingRightHandSide(Vec &input, Vec &output) = 0;
  
  virtual void initialize() = 0;
private:
};