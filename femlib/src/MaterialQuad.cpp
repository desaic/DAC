#include "MaterialQuad.hpp"
#include "Element.hpp"
#include "ElementMesh.hpp"
#include "femError.hpp"

#include "Quadrature.hpp"
#include "StrainEne.hpp"

#include <Eigen/Dense>

#include <iostream>

///@brief helper for computing stiffness contribution of one quadrature point
void stiffness(int qi, const MaterialQuad * mat, Element* ele, ElementMesh * mesh, Eigen::MatrixXd &K);

MaterialQuad::MaterialQuad(StrainEne * ene, Quadrature * _q ):q(_q)
{
  if(q==0){
    q = &(Quadrature::Gauss2);
  }
  e.resize(q->x.size(), ene);
}

MaterialQuad::MaterialQuad(const std::vector<StrainEne *> & ene, Quadrature * _q )
{
  if(q==0){
    q = &(Quadrature::Gauss2);
  }
  e=ene;
}

void MaterialQuad::init(ElementMesh * m)
{
}

double
MaterialQuad::getEnergy(int eidx, ElementMesh * mesh)
{
  double energy = 0;
  Element * ele = mesh->e[eidx];
  for(int ii = 0; ii<q->x.size();ii++){
    Eigen::Matrix3d F = ele->defGradCached(ii, mesh->X, mesh->x);
    if(mesh->dim==2){
      F(2,2)=1;
    }
    double det = F.determinant();
    if (det <= 0){
      fem_error = -1;
    }
    energy += ele->detJ[ii] * q->w[ii] * e[ii]->getEnergy(F);
  }

  if(mesh->dim == 2){
    energy *= mesh->depth[eidx];
  }
  return energy;
}

std::vector<Eigen::Vector3d>
MaterialQuad::getForce(int eidx, ElementMesh * mesh)
{
  Element * ele = mesh->e[eidx];
  int nV = ele->nV();
  std::vector<Eigen::Vector3d> f(nV, Eigen::Vector3d::Zero());
  std::vector<Eigen::Matrix3d> P(q->w.size());
  for(unsigned int ii = 0; ii<q->x.size(); ii++){
    Eigen::Matrix3d F = ele->defGradCached(ii, mesh->X, mesh->x);
    if(mesh->dim==2){
      F(2,2)=1;
    }
    if(F.determinant()<=0){
      fem_error = -1;
    }
    P[ii] = e[ii]->getPK1(F);
    if(mesh->dim == 2){
      P[ii](2,2) = 0;
    }
  }
  
  for(unsigned int jj = 0; jj<q->x.size(); jj++){
    for(int ii = 0; ii<nV; ii++){
      f[ii] -= ele->detJ[jj] * q->w[jj] * (P[jj] * ele->JdN[jj][ii]);
    }
  }

  if(mesh->dim==2){
    for(unsigned int ii =0 ; ii<f.size(); ii++){
      f[ii] *= mesh->depth[eidx];
    }
  }

  return f;
}

Eigen::MatrixXd
MaterialQuad::getStiffness(int eidx, ElementMesh * mesh)
{
  int dim = mesh->dim;

  Element * ele = mesh->e[eidx];
  int ndof = dim * ele->nV();
  Eigen::MatrixXd K = Eigen::MatrixXd::Zero(ndof, ndof);
  Eigen::MatrixXd Kq = Eigen::MatrixXd::Zero(ndof, ndof);

  for(unsigned int ii = 0; ii<q->x.size();ii++){
    Kq = Eigen::MatrixXd::Zero(ndof, ndof);
    stiffness(ii, this, ele, mesh, Kq);
    K += (ele->detJ[ii] * q->w[ii]) * Kq;
  }

  if(mesh->dim==2){
    K *= mesh->depth[eidx];
  }
  return K;
}

void
stiffness(int qi, const MaterialQuad * mat, Element* ele, ElementMesh * mesh,
          Eigen::MatrixXd & K)
{
  int dim = mesh->dim;
  int nV = ele->nV();

  Eigen::MatrixXd F = ele->defGradCached(qi, mesh->X,mesh->x);
  if(mesh->dim==2){
    F(2,2)=1;
  }

  for(int ii = 0;ii<nV;ii++){
    if(ii==0){
      mat->e[qi]->cacheF = false;
    }else{
      mat->e[qi]->cacheF = true;
    }
    std::vector<Eigen::Matrix3d> dP = mat->e[qi]->getdPdx(F, ele->JdN[qi][ii], dim);
    for(int vv = ii; vv<nV; vv++){
      for(int jj = 0;jj<dim;jj++){
        Eigen::Vector3d dfdxi = dP[jj] * ele->JdN[qi][vv];
        for(int kk = 0; kk<dim; kk++){
          K(dim*vv+kk, dim * ii + jj) = dfdxi[kk];
          K( dim * ii + jj, dim*vv+kk) = dfdxi[kk];
        }
      }
    }
  }
}

Eigen::MatrixXd
MaterialQuad::getMassMatrix(int eidx, ElementMesh * mesh)
{
  int dim = mesh->dim;
  Element * ele = mesh->e[eidx];
  int nV = ele->nV();
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(dim * nV, dim * nV);
  for(int qq = 0; qq<q->x.size(); qq++){
    Eigen::MatrixXd N = Eigen::MatrixXd::Zero(dim, dim*nV);
    std::vector<double> shape = ele->shapeFun(q->x[qq]);
    for(int ii = 0; ii<shape.size(); ii++){
      for(int jj = 0; jj<dim; jj++){
        N(jj, dim * ii + jj) = shape[ii];
      }
    }
    M += q->w[qq] * N.transpose() * N * ele->detJ[qq];
  }
  return M;
}

Eigen::Matrix3d MaterialQuad::getStress(int qidx, int eidx, ElementMesh * mesh)
{
  Element * ele = mesh->e[eidx];
  int nV = ele->nV();
  Eigen::Matrix3d P;
  Eigen::Matrix3d F = ele->defGradCached(qidx, mesh->X, mesh->x);
  if (mesh->dim == 2){
    F(2, 2) = 1;
  }
  if (F.determinant() <= 0){
    fem_error = -1;
  }
  P = e[qidx]->getPK1(F);
  if (mesh->dim == 2){
    P(2, 2) = 0;
  }
  return P;
}
