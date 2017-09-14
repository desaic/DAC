#ifndef STRAINCOROTLIN_HPP
#define STRAINCOROTLIN_HPP
#include "StrainLin.hpp"

#include <Eigen/Dense>

class StrainCorotLin :public StrainLin{
public:
  StrainCorotLin();
  virtual double getEnergy(const Eigen::Matrix3d & F);
  virtual Eigen::Matrix3d getPK1(const Eigen::Matrix3d & F);
  virtual std::vector<Eigen::Matrix3d>
    getdPdx(const Eigen::Matrix3d & F, const Eigen::Vector3d & dF, int dim=3);

};

#endif
