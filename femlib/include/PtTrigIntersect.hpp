#ifndef PT_TRIG_INTERSECT_HPP
#define PT_TRIG_INTERSECT_HPP
#include <Eigen/Dense>
struct InterResult
{
  bool intersect;
  double t;
  Eigen::Vector3d n;
  InterResult() :intersect(false), t(0){}
};

InterResult PtTrigIntersect(const Eigen::Vector3d & pt,
  const Eigen::Vector3d & dir,
  Eigen::Vector3d * trig);

///@param pt0 vertex in previous time step
///@param pt1 vertex in predicted time step
///@param trig0 an array of 3 vertices for a triangle in previous time step
///@param trig1 triangle in predicted time step
InterResult PtTrigIntersect(const Eigen::Vector3d & pt0,
  const Eigen::Vector3d & pt1,
  Eigen::Vector3d * trig0, Eigen::Vector3d * trig1);

#endif