#include "Element.hpp"
#include "femError.hpp"
#include <iostream>

std::vector<std::vector<Eigen::Vector3d> > Element::gradN;

double Element::getVol(const std::vector<Eigen::Vector3d> & X)
{
  return 0;
}

Eigen::Vector3d Element::getDisp(const Eigen::Vector3d & p,
                                 const std::vector<Eigen::Vector3d> & X,
                                 const std::vector<Eigen::Vector3d>x)
{
  std::vector<double> w = shapeFun(p);
  Eigen::Vector3d u(0,0,0);
  for(unsigned int ii = 0; ii<w.size(); ii++){
    u += w[ii]*(x[ii] - X[ii]);
  }
  return u;
}

Eigen::Matrix3d
Element::defGradCached(int qq,
                       const std::vector<Eigen::Vector3d> & X,
                       const std::vector<Eigen::Vector3d> & x) const
{
  Eigen::Matrix3d F = Eigen::Matrix3d::Zero();
  for(int ii = 0; ii<nV(); ii++){
    int vi = at(ii);
    //outer product
    F +=  x[vi] * gradN[qq][ii].transpose();
  }
  return F*Jinv[qq];
}

Eigen::Matrix3d
Element::getJ(int qq, const std::vector<Eigen::Vector3d> & X)
{
  Eigen::Matrix3d J = Eigen::Matrix3d::Zero();
  for(int ii = 0; ii<nV(); ii++){
    Eigen::Vector3d dN = gradN[qq][ii];
    //outer product
    J += X[at(ii)]* dN.transpose();
  }
//  std::cout<<J<<"\n";
  return J;
}

Element::Element(int _n):n(_n)
{}

Element::Element(const Element & e) : n(e.getNodeIndices())
{}
  

std::vector<std::array<int,2> >
Element::getEdges()
{
  std::vector<std::array<int,2> >edges;
  return edges;
}

void Element::resize(int size)
{
  n.resize(size);
}
