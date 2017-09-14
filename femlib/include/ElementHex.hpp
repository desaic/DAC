#ifndef ELEMENTHEX_HPP
#define ELEMENTHEX_HPP
#include "Element.hpp"

#include <Eigen/Dense>

extern const int facen[][3];
extern const int faces[][4];

class ElementHex:public Element
{
public:
  ElementHex();
  ElementHex(const ElementHex & e);

  int nF() const override{ return 6; }


  ///@brief natural coordinate for a point in reference space
  Eigen::Vector3d natCoord(const Eigen::Vector3d & p,
                    const std::vector<Eigen::Vector3d> & X);
  
  std::vector<double> shapeFun(const Eigen::Vector3d & p)const;

  Eigen::Vector3d shapeFunGrad(int ii, const Eigen::Vector3d &xx) const;

  std::vector<std::array<int,2> > getEdges();
  double getVol(const std::vector<Eigen::Vector3d> & X);

};
#endif
