#include "ElementHex.hpp"
#include "Quadrature.hpp"   //used to compute mass matrix
#include <iostream>

const int cubeEdges[12][2] =
{
  {0,1},
  {0,2},
  {1,3},
  {2,3},

  {4,5},
  {4,6},
  {5,7},
  {6,7},
  
  {0,4},
  {1,5},
  {2,6},
  {3,7}
};

static const int sw[8][3] =
{{-1,-1,-1},
 {-1,-1, 1}, 
 {-1, 1,-1},
 {-1, 1, 1},
 { 1,-1,-1},
 { 1,-1, 1},
 { 1, 1,-1},
 { 1, 1, 1}
};

const int faces[6][4]={
		{ 0, 1, 3, 2 },
		{ 4, 6, 7, 5 },
		{ 0, 4, 5, 1 },
		{ 2, 3, 7, 6 },
		{ 0, 2, 6, 4 },
		{ 1, 5, 7, 3 }
};

///@brief face normals
const int facen[6][3]={
	{-1, 0, 0},
	{ 1, 0, 0},
	{ 0,-1, 0},
	{ 0, 1, 0},
	{ 0, 0,-1},
	{ 0, 0, 1}
};

Eigen::Vector3d ElementHex::natCoord(const Eigen::Vector3d & p,
                                     const std::vector<Eigen::Vector3d> & X)
{
  Eigen::Vector3d n = p - X[at(0)];
  double size = (X[at(7)][0] - X[at(0)][0]);
  n = (2.0 / size) * n - Eigen::Vector3d(1,1,1);
  return n;
}

std::vector<double>
ElementHex::shapeFun(const Eigen::Vector3d & p)const
{
  std::vector<double> weights(8);
  for(int ii = 0;ii<nV();ii++){
    weights[ii] = (1.0f/8) * (1+sw[ii][0]*p[0])
      *(1+sw[ii][1]*p[1]) *(1+sw[ii][2]*p[2]) ;
  }
  return weights;
}

Eigen::Vector3d
ElementHex::shapeFunGrad(int ii, const Eigen::Vector3d & xx) const
{
  Eigen::Vector3d grad;
  grad[0] = 0.125*sw[ii][0] * (1 + sw[ii][1] * xx[1]) * (1 + sw[ii][2] * xx[2]);
  grad[1] = 0.125*sw[ii][1] * (1 + sw[ii][0] * xx[0]) * (1 + sw[ii][2] * xx[2]);
  grad[2] = 0.125*sw[ii][2] * (1 + sw[ii][0] * xx[0]) * (1 + sw[ii][1] * xx[1]);
  return grad;
}

std::vector<std::array<int,2> >
ElementHex::getEdges()
{
  int nEdge = 12;
  std::vector<std::array<int,2> >  edges(nEdge);
  for(int ii=0; ii<nEdge; ii++){
    edges[ii][0] = cubeEdges[ii][0];
    edges[ii][1] = cubeEdges[ii][1];
  }
  return edges;
}

double ElementHex::getVol(const std::vector<Eigen::Vector3d> & X)
{
  Eigen::Vector3d size = X[at(7)] - X[at(0)];
  double vol = size[0] * size[1] * size[2];
  return vol;
}

ElementHex::ElementHex():Element(8)
{}

ElementHex::ElementHex(const ElementHex & e) :Element(e)
{}
