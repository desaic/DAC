#include "StrainEneNeo.hpp"

StrainEneNeo::StrainEneNeo()
{
  param.resize(2);
  param[0] =1;
  param[1] =10;
}

double StrainEneNeo::getEnergy(const Eigen::Matrix3d &F)
{
  double I1 = (F.transpose()*F).trace();
  double JJ = std::log(F.determinant());
  double mu = param[0],lambda=param[1];
  double Psi = (mu/2) * (I1-3) - mu*JJ + (lambda/2)*JJ*JJ;
  return Psi;
}

Eigen::Matrix3d
StrainEneNeo::getPK1(const Eigen::Matrix3d & F)
{
  double JJ = std::log(F.determinant());
  Eigen::Matrix3d Finv = F.inverse().transpose();
  double mu = param[0],lambda=param[1];
  Eigen::Matrix3d PP = mu*(F-Finv) + lambda*JJ*Finv;
  return PP;
}

std::vector<Eigen::Matrix3d>
StrainEneNeo::getdPdx(const Eigen::Matrix3d & F,
                      const Eigen::Vector3d & dF,
                      int dim)
{
  if(!cacheF){
    double JJ = std::log(F.determinant());
    double mu = param[0];
    double lambda=param[1];
    c1 = mu-lambda * JJ;
    Eigen::Matrix3d FinvT = F.inverse().transpose();
    cacheFinvT = FinvT;
  }
  Eigen::Vector3d FinvTdF = cacheFinvT * dF;
  std::vector<Eigen::Matrix3d> dP(dim, Eigen::Matrix3d::Zero());
  for(int ii =0 ; ii<dim; ii++){
    dP[ii].row(ii) = param[0]*dF;
    dP[ii] += c1*FinvTdF * cacheFinvT.row(ii) + param[1] * FinvTdF[ii] * cacheFinvT;
  }
  return dP;
}
