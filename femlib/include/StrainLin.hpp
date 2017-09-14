#ifndef STRAIN_LIN_HPP
#define STRAIN_LIN_HPP
#include "StrainEne.hpp"
#include <Eigen/Dense>

class StrainLin :public StrainEne{
public:
	StrainLin();
  virtual double getEnergy(const Eigen::Matrix3d & F);
  virtual Eigen::Matrix3d getPK1(const Eigen::Matrix3d & F);
  virtual std::vector<Eigen::Matrix3d> getdPdx(const Eigen::Matrix3d & F,
                                  const Eigen::Vector3d & dF, int dim=3);

};

#endif
