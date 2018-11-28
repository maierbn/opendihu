#include "spatial_discretization/finite_element_method/solid_mechanics/solid_mechanics_utility.h"

#include <cmath>
#include <iostream>

namespace SpatialDiscretization
{

void checkSymmetry(double Cbar[3][3][3][3], std::string name)
{
  if (VLOG_IS_ON(1))
  {
    const double errorTolerance = 1e-14;
    for (int a=0; a<3; a++)
    {
      for (int b=0; b<3; b++)
      {
        for (int c=0; c<3; c++)
        {
          for (int d=0; d<3; d++)
          {
            if (fabs(Cbar[a][b][c][d] - Cbar[b][a][c][d]) > errorTolerance)
            {
              LOG(ERROR) << name << "[" <<a<< "][" <<b<< "][" << c<< "][" << d<< "] != " << name << "[" <<b<< "][" <<a<< "][" << c<< "][" << d<< "] (" <<Cbar[b][a][c][d]<< " != " <<Cbar[a][b][c][d]<< ") - minor symmetry violated" << std::endl;
            }
            if (fabs(Cbar[a][b][c][d] - Cbar[a][b][d][c]) > errorTolerance)
            {
              LOG(ERROR) << name << "[" <<a<< "][" <<b<< "][" << c<< "][" << d<< "] != " << name << "[" <<a<< "][" <<b<< "][" << d<< "][" << c<< "] (" <<Cbar[a][b][d][c]<< " != " <<Cbar[a][b][c][d]<< ") - minor symmetry violated" << std::endl;
            }
            if (fabs(Cbar[a][b][c][d] - Cbar[c][d][a][b]) > errorTolerance)
            {
              LOG(ERROR) << name << "[" <<a<< "][" <<b<< "][" << c<< "][" << d<< "] != " << name << "[" << c<< "][" << d<< "][" <<a<< "][" <<b<< "] (" <<Cbar[c][d][a][b]<< " != " <<Cbar[a][b][c][d]<< ") - major symmetry violated" << std::endl;
            }
          }
        }
      }
    }
  }
}

};  // namespace
