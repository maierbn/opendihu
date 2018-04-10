#include "spatial_discretization/finite_element_method/solid_mechanics/solid_mechanics_utility.h"

#include <cmath>
#include <iostream>

namespace SpatialDiscretization
{  

void checkSymmetry(double Cbar[3][3][3][3], std::string name)
{
  for (int a=0; a<3; a++)
  {
    for (int b=0; b<3; b++)
    {
      for (int c=0; c<3; c++)
      {
        for (int d=0; d<3; d++)
        {
          if (Cbar[a][b][c][d] != Cbar[b][a][c][d])
          {
            std::cout << name << "["<<a<<"]["<<b<<"]["<<c<<"]["<<d<<"] != " << name << "["<<b<<"]["<<a<<"]["<<c<<"]["<<d<<"] ("<<Cbar[b][a][c][d]<<" != "<<Cbar[a][b][c][d]<<") - minor symmetry violated"<<std::endl;
          }
          if (Cbar[a][b][c][d] != Cbar[a][b][d][c])
          {
            std::cout << name << "["<<a<<"]["<<b<<"]["<<c<<"]["<<d<<"] != " << name << "["<<a<<"]["<<b<<"]["<<d<<"]["<<c<<"] ("<<Cbar[a][b][d][c]<<" != "<<Cbar[a][b][c][d]<<") - minor symmetry violated"<<std::endl;
          }
          if (Cbar[a][b][c][d] != Cbar[c][d][a][b])
          {
            std::cout << name << "["<<a<<"]["<<b<<"]["<<c<<"]["<<d<<"] != " << name << "["<<c<<"]["<<d<<"]["<<a<<"]["<<b<<"] ("<<Cbar[c][d][a][b]<<" != "<<Cbar[a][b][c][d]<<") - major symmetry violated"<<std::endl;
          }
        }
      }
    }
  }
}
 
};  // namespace