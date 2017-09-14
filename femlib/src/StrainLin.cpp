#include "StrainLin.hpp"

StrainLin::StrainLin()
{
	param.resize(2);
	param[0] = 1;
	param[1] = 10;
}

double StrainLin::getEnergy(const Eigen::Matrix3d & F)
{
  Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d eps = 0.5*(F + F.transpose()) - I;
  double t = eps.trace();
  double Psi = param[0]*eps.squaredNorm() + 0.5f*param[1] * t*t;
	return Psi;
}

Eigen::Matrix3d
StrainLin::getPK1(const Eigen::Matrix3d & F)
{
  double mu = param[0];
  double lambda = param[1];
  Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d PP = mu*(F + F.transpose()-2*I) + lambda * (F.trace()-3) * I;
	return PP;
}

std::vector<Eigen::Matrix3d>
StrainLin::getdPdx(const Eigen::Matrix3d &F,
                   const Eigen::Vector3d &_dF, int dim)
{
  double mu = param[0];
  double lambda = param[1];
  Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
  std::vector<Eigen::Matrix3d> dP(dim);
  for(int ii =0 ; ii<dim; ii++){
    Eigen::Matrix3d dF = Eigen::Matrix3d::Zero();
    dF.row(ii) = _dF;
    dP[ii] = mu * (dF + dF.transpose()) + lambda * dF(ii,ii) * I;
  }

	return dP;
}
