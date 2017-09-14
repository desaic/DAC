#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include <array>
#include <vector>

#include <Eigen/Dense>

class Quadrature;

class Element
{
public:
  
  Element(int _n=0);
  Element(const Element & e);
  
  virtual ~Element(){}
  ///@brief number of vertices
  virtual int nV() const{return (int)n.size();}

  //@brief number of faces. Default returns 0 (not meaningful).
  virtual int nF() const{return 0;}

  int & operator[](int idx){return n[idx];}
  int operator[](int idx)const{return n[idx];}
  int at(int idx)const{return n[idx];}
  std::vector<int> verts()const { return n; }
  void resize(int size);

  ///@param qq quadrature index
  ///@param gradN precomputed shape function gradient
  ///@param X Reference vertex positions.
  ///@param x global array of Deformed vertex positions.
  Eigen::Matrix3d defGradCached(int qq,
    const std::vector<Eigen::Vector3d> & X,
    const std::vector<Eigen::Vector3d> & x) const;


  ///@param p natural coordinates.
  ///@return interpolation weight for each dof.
  virtual std::vector<double> shapeFun(const Eigen::Vector3d & p) const=0;

  ///@param ii index of basis function.
  ///@param xx point in natural coordinates.
  ///@param X global array of rest positions.
  virtual Eigen::Vector3d shapeFunGrad(int ii, const Eigen::Vector3d & xx) const=0;


  //pre-computed gradN for a unit-size element.
  //gradN[q][ii] is shape function gradient with respect to vertex ii and quadrature point q.
  //does not work for meshes with different element shape functions.
  static std::vector<std::vector<Eigen::Vector3d> > gradN;

  std::vector<std::vector<Eigen::Vector3d> > JdN;
  ///@brief inverse of material coordinate with respect to natural coordinate derivative.
  ///initialized by MaterialQuad::init
  std::vector<Eigen::Matrix3d> Jinv;
  
  ///@brief material space volume at each quadrature point.
  std::vector<double> detJ;

  ///@param p natural coordinate
  ///@param X rest pose.
  virtual Eigen::Matrix3d getJ(int qq, const std::vector<Eigen::Vector3d> & X);

  virtual std::vector<std::array<int,2> > getEdges();
  
  virtual double getVol(const std::vector<Eigen::Vector3d> & X);

  virtual Eigen::Vector3d getDisp(const Eigen::Vector3d & p,
                           const std::vector<Eigen::Vector3d> & X,
                           const std::vector<Eigen::Vector3d>x);

  const std::vector<int> & getNodeIndices()const{ return n; }
private:

  ///@brief nodal indices
  std::vector<int> n;

};
#endif
