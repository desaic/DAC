#include "ElementSpring.hpp" 

#include <iostream>

std::vector<double> ElementSpring::shapeFun(const Eigen::Vector3d & p)const
{
  std::vector<double> weights(nV(), 0);
  weights[0] = (1.0f/2) * (  p[0] + 1);
  weights[1] = (1.0f/2) * (- p[0] + 1);
  return weights;
}

Eigen::Vector3d
ElementSpring::shapeFunGrad(int ii, const Eigen::Vector3d & xx,
                                 const std::vector<Eigen::Vector3d> & X) const
{
  //rest length
  double len=(X[at(1)] - X[at(0)]).norm();
//  std::cout<<size[0]<<"\n";
  Eigen::Vector3d grad;
  len = 1.0/len;
  grad[0] = 1;
  grad[1] = 0;
  grad[2] = 0;
  return grad;
}

std::vector<std::array<int,2> >
ElementSpring::getEdges()
{
  int nEdge = 1;
  std::vector<std::array<int,2> >  edges(nEdge);
  edges[0][0] = 0;
  edges[0][1] = 1;
  return edges;
}

double ElementSpring::getVol(const std::vector<Eigen::Vector3d> & X)
{
  Eigen::Vector3d size = X[at(3)] - X[at(0)];
  double vol = size.norm();
  return vol;
}

ElementSpring::ElementSpring():Element(2)
{}

ElementSpring::ElementSpring(const ElementSpring & e) :Element(e)
{}

