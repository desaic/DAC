#include "PtTrigIntersect.hpp"
#include "CTCD.hpp"

InterResult PtTrigIntersect(const Eigen::Vector3d & pt,
  const Eigen::Vector3d & dir,
  Eigen::Vector3d * trig)
{
  InterResult result;
  Eigen::Vector3d diff = pt - trig[0];
  Eigen::Vector3d edge1 = trig[1] - trig[0];
  Eigen::Vector3d edge2 = trig[2] - trig[0];
  Eigen::Vector3d normal = edge1.cross(edge2);
  Eigen::Vector3d direction = normal;  //using normal as direction is better so far
  //using linear velocity of vertices will miss collisions a lot of times.
  //ideally should use vertex cube intersection code.
  // Solve Q + t*D = b1*E1 + b2*E2 (Q = diff, D = line direction,
  // E1 = edge1, E2 = edge2, N = Cross(E1,E2)) by
  //   |Dot(D,N)|*b1 = sign(Dot(D,N))*Dot(D,Cross(Q,E2))
  //   |Dot(D,N)|*b2 = sign(Dot(D,N))*Dot(D,Cross(E1,Q))
  //   |Dot(D,N)|*t = -sign(Dot(D,N))*Dot(Q,N)
  double DdN = direction.dot(normal);

  double DdQxE2 = direction.dot(diff.cross(edge2));
  if (DdQxE2 >= 0)
  {
    double DdE1xQ = direction.dot(edge1.cross(diff));
    if (DdE1xQ >= 0)
    {
      if (DdQxE2 + DdE1xQ <= DdN)
      {
        // Line intersects triangle.
        double QdN = -diff.dot(normal);
        double inv = 1 / DdN;

        result.intersect = true;
        result.t = QdN*inv;
        result.n = normal.normalized();
        return result;
      }
      // else: b1+b2 > 1, no intersection
    }
    // else: b2 < 0, no intersection
  }
  // else: b1 < 0, no intersection

  result.intersect = false;
  return result;
}

InterResult PtTrigIntersect(const Eigen::Vector3d & pt0,
  const Eigen::Vector3d & pt1,
  Eigen::Vector3d * trig0, Eigen::Vector3d * trig1)
{
  Eigen::Vector3d edge1 = trig1[1] - trig1[0];
  Eigen::Vector3d edge2 = trig1[2] - trig1[0];
  Eigen::Vector3d normal = edge1.cross(edge2);
  // TODO:
  // 1. didn't have time to fix depenedencies in compilation so you'll need to compile and test
  // 2. Note that all values are passed in by Eigen::Vector3d while you're using float--might want to template out float vs double everywhere generically either way resovle here.
  // 3. See below for todos for setting result.t and result.n
  InterResult result;
  double eta = 1e-4; // the distance from coplanarity when we declare collision
  double t = 0; // earliest detected time of collision on interval
  result.intersect = CTCD::vertexFaceCTCD(pt0, trig0[0], trig0[1], trig0[2], pt1, trig1[0], trig1[1], trig1[2], eta, t);
  if (result.intersect){
    result.t = t; // TODO: warning! I'm overloading your t with 1st time of collision
    result.n = normal.normalized(); //TODO: as this is the normal you are using for contact resolution this should be the outward pointing face normal--here I'm using end of postion normal alternately you could use start normal as above or normal at time t.
  }
  return result;
}
