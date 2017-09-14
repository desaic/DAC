#include "StrainEneANeo.hpp"

StrainEneANeo::StrainEneANeo()
{
  param.resize(6);
  param[0] =1;
  param[1] =10;
  param[2] = 1;
  param[3] = 1;
  param[4] = 1;
  param[5] = 1;
  double s2 = 0.5*std::sqrt(2);
  int nSpring = 4;
  //x y z axis
  V.resize(nSpring, Eigen::Vector3d::Zero());
  for (unsigned int ii = 0; ii<V.size(); ii++){
    //V[ii][ii] = 1;
    //V[3+ii] = s2 * Eigen::Vector3d::Ones();
    //V[3 + ii][ii] = 0;
  }
  V[0] = Eigen::Vector3d(s2, s2, 0);
  V[1] = Eigen::Vector3d(-s2, s2, 0);
  V[2] = Eigen::Vector3d(0, s2, s2);
  V[3] = Eigen::Vector3d(0, s2, -s2);
}

double StrainEneANeo::getEnergy(const Eigen::Matrix3d &F)
{
  double I1 = (F.transpose()*F).trace();
  double JJ = std::log(F.determinant());
  double mu = param[0], lambda = param[1];
  double Psi = (mu / 2) * (I1 - 3) - mu*JJ + (lambda / 2)*JJ*JJ;

  for (unsigned int ii = 0; ii<V.size(); ii++){
    Eigen::Vector3d spring = F*V[ii];
    double len = (spring.norm() - 1);
    Psi += 0.5 * param[2 + ii] * len * len;
  }

  //Eigen::Vector3d spring = F*V[1];
  //double J1 = std::log(spring.squaredNorm());
  //Psi += 0.5 * param[3] * J1 * J1;

  return Psi;
}

Eigen::Matrix3d
StrainEneANeo::getPK1(const Eigen::Matrix3d & F)
{
  double JJ = std::log(F.determinant());
  Eigen::Matrix3d Finv = F.inverse().transpose();
  double mu = param[0], lambda = param[1];
  Eigen::Matrix3d PP = mu*(F - Finv) + lambda*JJ*Finv;

  for (unsigned int ii = 0; ii<V.size(); ii++){
    Eigen::Vector3d spring = F*V[ii];
    float len = spring.norm();
    float coef = param[ii + 2] * (len - 1) / len;
    PP += coef * spring * V[ii].transpose();
  }

  //Eigen::Vector3d spring = F*V[1];
  //double l = spring.squaredNorm();
  //double J1 = std::log(spring.squaredNorm());
  //Eigen::Matrix3d P1 = (2 *J1 / l) * F * (V[1] * V[1].transpose());
  //PP += param[3] * P1;

  return PP;
}

std::vector<Eigen::Matrix3d>
StrainEneANeo::getdPdx(const Eigen::Matrix3d & F,
                      const Eigen::Vector3d & dF,
                      int dim)
{
  if (!cacheF){
    double JJ = std::log(F.determinant());
    double mu = param[0];
    double lambda = param[1];
    c1 = mu - lambda * JJ;
    Eigen::Matrix3d FinvT = F.inverse().transpose();
    cacheFinvT = FinvT;
  }
  Eigen::Vector3d FinvTdF = cacheFinvT * dF;
  std::vector<Eigen::Matrix3d> dP(dim, Eigen::Matrix3d::Zero());
  for (int ii = 0; ii<dim; ii++){
    dP[ii].row(ii) = param[0] * dF;
    dP[ii] += c1*FinvTdF * cacheFinvT.row(ii) + param[1] * FinvTdF[ii] * cacheFinvT;
  }

  for (unsigned int ii = 0; ii<V.size(); ii++){
    Eigen::Vector3d spring = F*V[ii];
    float len = spring.norm();
    float coef1 = param[ii + 2] * (len - 1) / (len*len*len);
    float coef2 = param[ii + 2] / (len*len*len);
    Eigen::Matrix3d dPdSpring = coef1 * (len*len) * Eigen::Matrix3d::Identity()
      + coef2*spring*spring.transpose();
    for (int di = 0; di < dim; di++){
      Eigen::Vector3d dSpring = Eigen::Vector3d::Zero();
      dSpring[di] = dF.transpose()*V[ii];
      //derivative of spring outer V when spring changed dSpring.
      Eigen::Matrix3d dSpringV = dSpring * V[ii].transpose();
      dP[di] += dPdSpring * dSpringV;
    }
  }

  //Eigen::Vector3d spring = F*V[1];
  //double l = spring.squaredNorm();
  //double J1 = std::log(spring.squaredNorm());
  //Eigen::Matrix3d VVT = V[1] * V[1].transpose();
  //for (int di = 0; di < dim; di++){
  //  Eigen::Matrix3d dF1 = Eigen::Matrix3d::Zero();
  //  dF1.row(di) = dF;
  //  double dl = 2*(F*VVT).row(di) * dF;
  //  Eigen::Matrix3d dP1 = (2*dl * (1- J1)/ (l*l) ) * F * VVT  + (2*J1/l) * dF1 * VVT;
  //  dP[di] += param[3] * dP1;
  //}

  return dP;
}
