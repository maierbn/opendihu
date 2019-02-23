
#ifndef TRANSFER_H
#define TRANSFER_H
namespace Transfer
{
class Transfer {
    
public:
    
//!Restriction step of the multigrid method
void restriction(std::vector<double> *valuesFine, std::vector<double> *valuesCoarse);

//!prolongation step of the multigrid method
void prolongation(std::vector<double> *valuesCoarse, std::vector<double> *valuesFine);

};

}
#include "transfer.tpp"
#endif //TRANSFER_H
