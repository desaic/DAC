#ifndef BOX_HPP
#define BOX_HPP

#include <Eigen/Dense>
#include <vector>
struct Box{
  Eigen::Vector3d mn, mx;
  Box() :mn(0, 0, 0), mx(0, 0, 0){}
  //points samples on each edge.
  std::vector<Eigen::Vector3d> edgeVerts;
  void computeEdgeVerts(){
    //num vert on each edge.
    int N = 10;
    for (int i = 0; i < N; i++){
      //linear interpolation weight
      double w = (i + 1) / (N + 1);
      Eigen::Vector3d v0 = mn, v1 = mn;
      v0[1] = mx[1];
      v1[1] = mx[1];
      v1[2] = mx[2];
      edgeVerts.push_back((1 - w)*v0 + w*v1);
      v0[0] = mx[0];
      v1[0] = mx[0];
      edgeVerts.push_back((1 - w)*v0 + w*v1);
    }
  }
};
#endif
