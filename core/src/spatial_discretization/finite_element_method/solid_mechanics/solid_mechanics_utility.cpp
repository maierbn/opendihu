#include "spatial_discretization/finite_element_method/solid_mechanics/solid_mechanics_utility.h"

#include <cmath>
#include <iostream>

namespace SpatialDiscretization
{  

void checkSymmetry(double Cbar[3][3][3][3], std::string name)
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
            LOG(ERROR) << name << "["<<a<<"]["<<b<<"]["<<c<<"]["<<d<<"] != " << name << "["<<b<<"]["<<a<<"]["<<c<<"]["<<d<<"] ("<<Cbar[b][a][c][d]<<" != "<<Cbar[a][b][c][d]<<") - minor symmetry violated"<<std::endl;
          }
          if (fabs(Cbar[a][b][c][d] - Cbar[a][b][d][c]) > errorTolerance)
          {
            LOG(ERROR) << name << "["<<a<<"]["<<b<<"]["<<c<<"]["<<d<<"] != " << name << "["<<a<<"]["<<b<<"]["<<d<<"]["<<c<<"] ("<<Cbar[a][b][d][c]<<" != "<<Cbar[a][b][c][d]<<") - minor symmetry violated"<<std::endl;
          }
          if (fabs(Cbar[a][b][c][d] - Cbar[c][d][a][b]) > errorTolerance)
          {
            LOG(ERROR) << name << "["<<a<<"]["<<b<<"]["<<c<<"]["<<d<<"] != " << name << "["<<c<<"]["<<d<<"]["<<a<<"]["<<b<<"] ("<<Cbar[c][d][a][b]<<" != "<<Cbar[a][b][c][d]<<") - major symmetry violated"<<std::endl;
          }
        }
      }
    }
  }
}

void checkSymmetry(const std::array<Vec3,3> &rightCauchyGreen, std::string name)
{
  bool isSymmetric = true;
  const double errorTolerance = 1e-14;
  for (int a=0; a<3; a++)
  {
    for (int b=0; b<3; b++)
    {
      if (fabs(rightCauchyGreen[a][b] - rightCauchyGreen[b][a]) > errorTolerance)
      {
        LOG(ERROR) << name << "["<<a<<"]["<<b<<"] != " << name << "["<<b<<"]["<<a<<"] ("<<rightCauchyGreen[a][b]<<" != "<<rightCauchyGreen[b][a]<<") - symmetry violated"<<std::endl;
        isSymmetric = false;
      }
    }
  }
  if (isSymmetric)
    VLOG(2) << name << " is symmetric";
}

void checkInverseIsCorrect(const std::array<Vec3,3> &rightCauchyGreen, const std::array<Vec3,3> &inverseRightCauchyGreen, std::string name)
{
  bool inverseCorrect = true;
  double errorTolerance = 1e-14;
  for (int a=0; a<3; a++)
  {
    for (int b=0; b<3; b++)
    {
      double matrixProduct = 0.0;
      for (int k=0; k<3; k++)
      {
        matrixProduct += rightCauchyGreen[k][a] * inverseRightCauchyGreen[b][k];
      }
      double delta_ab = (a == b? 1.0 : 0.0);
      
      if (fabs(delta_ab - matrixProduct) > errorTolerance)
      {
        LOG(ERROR) << name << " or " << name <<"^{-1} is wrong: " << matrixProduct << " should be " << delta_ab;
        inverseCorrect = false;
      }
    }
  }
  
  if (inverseCorrect)
  {
    VLOG(2) << name << " corresponds to " << name <<"^{-1}";
  }
  else
  {
    LOG(DEBUG) << name << ": " << rightCauchyGreen << std::endl << ", inv: " << inverseRightCauchyGreen << std::endl;
  }
}

 
};  // namespace