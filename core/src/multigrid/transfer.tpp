#include "multigrid/transfer.h"

#include "utility/python_utility.h"

namespace Transfer
{

//!Restriction step of the multigrid method
void restriction(std::vector<double> *valuesFine, std::vector<double> *valuesCoarse){
	valuesFine = valuesCoarse;
}
//!prolongation step of the multigrid method
void prolongation(std::vector<double> *valuesCoarse, std::vector<double> *valuesFine){
	valuesCoarse = valuesFine;
}
};
