#include "WalkerObjective.hpp"
#include "MaterialQuad.hpp"

void WalkerObjective::init(const Eigen::VectorXd & x0)
{
  materials.resize(em->m.size());
  for (unsigned int ii = 0; ii<em->m.size(); ii++){
    materials[ii] = (MaterialQuad*)em->m[ii];
  }
}

///@brief set simulation to time 0. Restore external forces, x, v, etc.
void WalkerObjective::resetSim()
{
  world->resetSim();
}

///@override
///@param x. 6 parameters.
//x, y z size for body and leg 
void WalkerObjective::setParam(const Eigen::VectorXd & x0)
{
  param = x0;
  for (int i = 0; i < x0.rows(); i++){
    std::cout << x0[i] << " ";
  }
  std::cout << "\n";

  walkerParam.bodySize[0] = x0[0];
  walkerParam.bodySize[1] = x0[1];
  walkerParam.bodySize[2] = x0[2];

  walkerParam.legSize[0] = x0[3];
  walkerParam.legSize[1] = x0[4];
  walkerParam.legSize[2] = x0[5];

  em->vertIdx.clear();
  makeWalker(em, &walkerParam);
  em->initElements(materials[0]->q);
  em->Kpattern.resize(0, 0); 
  em->sblocks.clear();
  Eigen::SparseMatrix<double> K = em->getStiffnessSparse();
  world->solvers[0]->init(K);
}

///@override
double WalkerObjective::f()
{
  resetSim();

  std::string filename = "opt_mesh" + std::to_string(nEval) + ".txt";
  std::cout << "save mesh file " << filename << "\n";
  std::ofstream out(filename);
  em->saveMesh(out);
  out.close();

  world->loop();
  //minimize negative speed.
  double val = 1-world->overall_speed;
  if (val != val){
    val = 100;
  }

  std::cout << "obj val "<< nEval<<" "<< val << "\n";

  nEval++;
  return val;
}

Eigen::VectorXd WalkerObjective::df()
{
  Eigen::VectorXd grad(param.size());
  Eigen::VectorXd x0 = param;
  for (size_t i = 0; i < param.size(); i++){
    x0[i] += dParam;
    setParam(x0);
    double fp = f();

    x0[i] -= 2 * dParam;
    setParam(x0);
    double fm = f();

    x0[i] += dParam;

    grad[i] = (fp - fm) / (2 * dParam);
  }
  param = x0;
  return grad;
}
