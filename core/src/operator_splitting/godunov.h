#pragma once

#include "operator_splitting/operator_splitting.h"
#include "control/runnable.h"

namespace OperatorSplitting
{

template<typename Runnable1, typename Runnable2>
class Godunov : public OperatorSplitting, Runnable
{
public:
private:
};

}  // namespace