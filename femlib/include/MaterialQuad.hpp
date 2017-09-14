#ifndef MATERIALQUAD_HPP
#define MATERIALQUAD_HPP

#include "Material.hpp"
#include <Eigen/Dense>

class StrainEne;
class Quadrature;

///@brief abstract material model that uses quadrature
class MaterialQuad:public Material
{
public:
  MaterialQuad(StrainEne * ene = 0, Quadrature * _q = 0);

  MaterialQuad(const std::vector<StrainEne *> &ene, Quadrature * _q = 0);
  
  ///@brief precompute gradN, Jinv, detJ.
  void init(ElementMesh * mesh);

  double getEnergy(int eidx, ElementMesh * mesh);
  std::vector<Eigen::Vector3d> getForce(int eidx, ElementMesh * mesh);
  
  Eigen::MatrixXd getStiffness(int eidx, ElementMesh * mesh);

  virtual Eigen::MatrixXd getMassMatrix(int eidx, ElementMesh * mesh);
  
  virtual Eigen::Matrix3d getStress(int qidx, int eidx, ElementMesh * mesh);
  
  std::vector<StrainEne*> e;
  const Quadrature * q;

};

#endif
