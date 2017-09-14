#ifndef STRAINENE_HPP
#define STRAINENE_HPP

#include <Eigen/Dense>

#include <vector>

///@brief abstract class for strain energy functions, i.e. function of deformation gradient
///also computes derivatives
class StrainEne{
public:
  StrainEne():cacheF(0){}
  virtual double getEnergy(const Eigen::Matrix3d & F)=0;
  virtual Eigen::Matrix3d getPK1(const Eigen::Matrix3d & F)=0;
  virtual std::vector<Eigen::Matrix3d>
  getdPdx(const Eigen::Matrix3d & F, const Eigen::Vector3d & dF, int dim=3)=0;
  virtual ~StrainEne();
  
  std::vector<double> param;

  enum MaterialModel{
    LIN, COROT, NEO, ANEO
  };

  ///@brief cache quantities determined by F.
  bool cacheF;
};

std::vector<StrainEne* > loadMaterials(std::istream & in,
  std::vector<double> & densities);

#endif
