#ifndef MATERIAL_HPP
#define MATERIAL_HPP
#include <vector>
#include <Eigen/Dense>

class Element;
class ElementMesh;

class Material
{
public:
  Material();
  ///@param eidx Element index
  virtual double getEnergy(int eidx, ElementMesh * mesh)=0;
  virtual std::vector<Eigen::Vector3d> getForce(int eidx, ElementMesh * mesh)=0;
  virtual Eigen::MatrixXd getStiffness(int eidx, ElementMesh * mesh);
  virtual Eigen::MatrixXd getMassMatrix(int eidx, ElementMesh * mesh);
  ///\brief get stress at quadrature points.
  /// @param qidx quadrature index.
  /// @param eidx Element index.
  virtual Eigen::Matrix3d getStress(int qidx, int eidx, ElementMesh * mesh);
  virtual ~Material();

  virtual void init(ElementMesh * mesh){}
  
  std::vector<float> param;
  float density;
};

#endif
