#include "ElementQuad.hpp"
#include <iostream>

const int quadEdges[4][2] =
{ {0,1},
  {1,3},
  {3,2},
  {2,0}};

const int sw[4][2] =
{{ -1,-1},
 { -1, 1},
 {  1,-1},
 {  1, 1}};

Eigen::Vector3d
ElementQuad::natCoord(const Eigen::Vector3d &p,
                      const std::vector<Eigen::Vector3d> &X)
{
  Eigen::Vector3d n = p - X[at(0)];
  float size = (X[at(3)][0] - X[at(0)][0]);
  n = (2.0 / size) * n - Eigen::Vector3d(1,1,1);
  return n;
}

std::vector<double>
ElementQuad::shapeFun(const Eigen::Vector3d & p)const
{
  std::vector<double> weights(nV(), 0);
  for(int ii = 0;ii<nV();ii++){
    weights[ii] = 0.25 * (1+sw[ii][0]*p[0])*(1+sw[ii][1]*p[1]);
  }
  return weights;
}

Eigen::Vector3d
ElementQuad::shapeFunGrad(int ii, const Eigen::Vector3d  & p) const
{
  Eigen::Vector3d grad;
  grad[0] = 0.25 * sw[ii][0] * (1 + sw[ii][1] * p[1]);
  grad[1] = 0.25 * sw[ii][1] * (1 + sw[ii][0] * p[0]);
  grad[2] = 0;
  return grad;
}

std::vector<std::array<int,2> >
ElementQuad::getEdges()
{
  int nEdge = 4;
  std::vector<std::array<int,2> >  edges(nEdge);
  for(int ii=0; ii<nEdge; ii++){
    edges[ii][0] = quadEdges[ii][0];
    edges[ii][1] = quadEdges[ii][1];
  }
  return edges;
}

double ElementQuad::getVol(const std::vector<Eigen::Vector3d > & X)
{
  Eigen::Vector3d size = X[at(3)] - X[at(0)];
  double vol = size[0] * size[1];
  return vol;
}

ElementQuad::ElementQuad():Element(4)
{}

ElementQuad::ElementQuad(const ElementQuad & e) :Element(e)
{}
