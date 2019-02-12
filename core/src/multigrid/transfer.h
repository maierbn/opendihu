
#ifndef TRANSFER_H
#define TRANSFER_H
namespace Transfer
{
class Transfer {
    
public:
    
//!Restriction step of the multigrid method
void restriction(std::vector<double> *values);

//!prolongation step of the multigrid method
void prolongation(std::vector<double> *values);

};

}
#include "transfer.tpp"
#endif //TRANSFER_H
