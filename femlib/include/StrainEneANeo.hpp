#ifndef STRAINENEANEO_HPP
#define STRAINENEANEO_HPP
#include "StrainEne.hpp"
#include <Eigen/Dense>
//anisotropic neohookean

class StrainEneANeo:public StrainEne
{
public:
  StrainEneANeo();
  double getEnergy(const Eigen::Matrix3d & F);
  Eigen::Matrix3d getPK1(const Eigen::Matrix3d & F);
  std::vector<Eigen::Matrix3d>
  getdPdx(const Eigen::Matrix3d &F, const Eigen::Vector3d &dF, int dim=3);

  Eigen::Matrix3d cacheFinvT;
  double c1;
  std::vector<Eigen::Vector3d> V;
};
#endif
