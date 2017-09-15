#include "CatapultObjective.hpp"
#include "WorldC.hpp"

CatapultObjective::CatapultObjective()
  :dParam(1e-3)
{}

void CatapultObjective::init(const Eigen::VectorXd & x0)
{
  
  conf->getFloat("beamwidth", opts.beamwidth);
  conf->getFloat("leftmargin", opts.leftmargin);
  segmentJumper3d( *(world->em_[0]), opts);
  X0 = world->em_[0]->X;
  param = x0;
  //set initial point to streth in config file
  std::vector<float> stretch = conf->getFloatVector("stretch");
  if (stretch.size() == 3){
    opts.offset = Eigen::Vector3d(stretch[0], stretch[1], stretch[2]);
    param[0] = stretch[0];
    param[1] = stretch[1];
  }
}

///@override
///@param x. Two parameters. Width and height.
void CatapultObjective::setParam(const Eigen::VectorXd & x0)
{
  opts.offset[0] = x0[0];
  opts.offset[1] = x0[1];
  world->em_[0]->X = X0;
  stretchJumper3d(*(world->em_[0]), opts);
  world->resetSim();
  world->loop();
}

///@override
double CatapultObjective::f()
{
  Eigen::Vector3d c = world->em_[1]->centerOfMass();
  Eigen::Vector3d p = world->em_[1]->linearMomentum();
  c[0] -= (world->catapultBox.mx[0] + world->wallDist);
  c[1] -= world->targetHeight;
  double w_p = 1e-4;
  double obj = c.squaredNorm() + w_p * p[0];
  return obj;
}
Eigen::VectorXd 
CatapultObjective::df()
{
  Eigen::VectorXd grad(param.size());
  Eigen::VectorXd x0 = param;
  for (size_t i = 0; i < param.size(); i++){
    x0[i] += dParam;
    setParam(x0);
    double fp = f();

    x0[i] -= 2*dParam;
    setParam(x0);
    double fm = f();

    x0[i] += dParam;

    grad[i] = (fp - fm) / (2*dParam);
  }
  param = x0;
  return grad;
}
