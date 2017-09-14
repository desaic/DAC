#include "StrainCorotLin.hpp"
#include <Eigen/Dense>
#include <vector>

#include <iostream>

std::vector<Eigen::Matrix3d> initEp3();
//3D Levi-Civita symbol
std::vector<Eigen::Matrix3d> Ep3=initEp3();

std::vector<Eigen::Matrix3d> initEp3()
{
  std::vector<Eigen::Matrix3d> Ep3(3, Eigen::Matrix3d::Zero());
  Ep3[0](1, 2) = 1;
  Ep3[0](2, 1) = -1;

  Ep3[1](0, 2) = -1;
  Ep3[1](2, 0) = 1;

  Ep3[2](0, 1) = 1;
  Ep3[2](1, 0) = -1;

  return Ep3;
}

double dot(Eigen::Matrix3d m1, Eigen::Matrix3d m2){
  double prod = 0;
  for (int ii = 0; ii < 3; ii++){
    for (int jj = 0; jj < 3; jj++){
      prod += m1(ii, jj) * m2(ii, jj);
    }
  }
  return prod;
}

StrainCorotLin::StrainCorotLin()
{
  param.resize(2);
  param[0] = 1;
  param[1] = 10;
}

double StrainCorotLin::getEnergy(const Eigen::Matrix3d &F)
{
  Eigen::JacobiSVD<Eigen::Matrix3d> svd(F);
  Eigen::Vector3d Sigma = svd.singularValues() - Eigen::Vector3d(1,1,1);
  double mu = param[0];
  double lambda = param[1];
  Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
  double t = Sigma[0] + Sigma[1] + Sigma[2];
  double Psi = mu*(Sigma[0] * Sigma[0] + Sigma[1] * Sigma[1] + Sigma[2] * Sigma[2]) + 0.5f * lambda * t * t;
  return Psi;
}

Eigen::Matrix3d StrainCorotLin::getPK1(const Eigen::Matrix3d & F)
{
  Eigen::JacobiSVD<Eigen::Matrix3d> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Eigen::Matrix3d U = svd.matrixU();
  Eigen::Matrix3d V = svd.matrixV();
  Eigen::Matrix3d R = U * V.transpose();
  double mu = param[0];
  double lambda = param[1];
  Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d P = 2*mu*(F-R) + lambda * (R.transpose() * F - I).trace() * R;

  return P;
}

Eigen::Matrix3d crossProdMat(const Eigen::Vector3d & v)
{
  Eigen::Matrix3d A = Eigen::Matrix3d::Zero();
  A(0, 1) = -v[2];
  A(0, 2) = v[1];
  A(1, 0) = v[2];
  A(1, 2) = -v[0];
  A(2, 0) = -v[1];
  A(2, 1) = v[0];
  return A;
}

std::vector<Eigen::Matrix3d>
StrainCorotLin::getdPdx(const Eigen::Matrix3d & F, const Eigen::Vector3d &_dF
                        , int dim)
{
  Eigen::JacobiSVD<Eigen::Matrix3d> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Eigen::Matrix3d U = svd.matrixU();
  Eigen::Matrix3d V = svd.matrixV();
  Eigen::Matrix3d R = U * V.transpose();
  Eigen::Matrix3d Sigma = svd.singularValues().asDiagonal();
  Eigen::Matrix3d S = V*Sigma*V.transpose();
  Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d SI = (S.trace()*I - S).inverse();

  //debug dR computation
  //double h = 0.01;
  //Eigen::Matrix3d F1 = F + h*dF;
  //svd.compute(F1, Eigen::ComputeFullU | Eigen::ComputeFullV);
  //U = svd.matrixU();
  //V = svd.matrixV();
  //Eigen::Matrix3d R1 = U * V.transpose();
  //R1 -= R;
  //R1 = (1/h)*R1;
  //std::cout << dR << "\n" << R1 << "\n\n";

  double mu = param[0];
  double lambda = param[1];

  std::vector<Eigen::Matrix3d> dP(dim);

  for(int ii =0 ; ii<dim; ii++){
    Eigen::Matrix3d dF=Eigen::Matrix3d::Zero();
    dF.row(ii) = _dF;
    Eigen::Matrix3d W = R.transpose()*dF;
    Eigen::Vector3d w;
    w[0] = W(1, 2) - W(2, 1);
    w[1] = W(2, 0) - W(0, 2);
    w[2] = W(0, 1) - W(1, 0);
    Eigen::Matrix3d dR = -R*crossProdMat(SI*w);
    dP[ii] = 2 * mu*dF + lambda*W.trace()*R +(lambda*(S - I).trace() - 2 * mu)*dR;
  }
  return dP;
}
