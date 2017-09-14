#ifndef ELEMENTQUAD_HPP
#define ELEMENTQUAD_HPP
#include "Element.hpp"

#include <Eigen/Dense>

class ElementQuad:public Element
{
public:
  ElementQuad();
  ElementQuad(const ElementQuad & e);

  int nF() const override{ return 1; }

  std::vector<std::array<int,2> > getEdges();

  ///@brief natural coordinate for a point in reference space
  Eigen::Vector3d natCoord(const Eigen::Vector3d & p,
                           const std::vector<Eigen::Vector3d> & X);

  ///@param p parameter in natural coodinates
  std::vector<double> shapeFun(const Eigen::Vector3d &p)const;

  ///@brief derivative with respect to material coodinate.
  Eigen::Vector3d shapeFunGrad(int ii, const Eigen::Vector3d &xx) const;

  double getVol(const std::vector<Eigen::Vector3d> &X);

};
#endif
