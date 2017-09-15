#include "Tests.hpp"
#include <Eigen/Dense>
#include "CTCD.hpp"
#include "WorldC.hpp"
#include "PtTrigIntersect.hpp"

void testCTCD(){
  Eigen::Vector3d x0, x1, trig0[3], trig1[3];
  x0 = Eigen::Vector3d(1.0, 1.0, 1.0);
  x1 = Eigen::Vector3d(0.9, 0.9, 0.9);
  trig0[0] = Eigen::Vector3d(1, 0, 0);
  trig0[1] = Eigen::Vector3d(0, 1, 0);
  trig0[2] = Eigen::Vector3d(0, 0, 1);
  
  trig1[0] = Eigen::Vector3d(3, 0, 0);
  trig1[1] = Eigen::Vector3d(0, 3, 0);
  trig1[2] = Eigen::Vector3d(0, 0, 3);

  InterResult result = PtTrigIntersect(x0, x1, trig0, trig1);

  std::cout << result.intersect<<"\n";
  std::cout << result.t<<" "<< (x0 + result.t * (x1-x0)) <<"\n";

  std::cout << result.n;
}
