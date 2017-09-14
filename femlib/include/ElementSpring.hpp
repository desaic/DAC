#ifndef ELEMENTSPRING_HPP
#define ELEMENTSPRING_HPP

#include "Element.hpp"

class ElementSpring: public Element{
public:
  ElementSpring();
  ElementSpring(const ElementSpring & e);

  int nF() const override{ return 0; }

  std::vector<std::array<int,2> > getEdges();

  std::vector<double> shapeFun(const Eigen::Vector3d & p)const;

  Eigen::Vector3d shapeFunGrad(int ii, const Eigen::Vector3d & xx,
    const std::vector<Eigen::Vector3d> & X) const;

  double getVol(const std::vector<Eigen::Vector3d> & X);
};

#endif
