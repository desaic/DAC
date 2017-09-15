#include "ArrayUtil.hpp"
#include "EigenUtil.hpp"
#include "ElementQuad.hpp"
#include "FileUtil.hpp"
#include "JumperObjective.hpp"
#include "LinPardiso.hpp"
#include "MaterialQuad.hpp"
#include "Quadrature.hpp"
#include "StepperNewmark.hpp"
#include "StrainEne.hpp"
#include "StretchMesh.hpp"
#include "World.hpp"
#include "Stepper.hpp"
#include <algorithm>
#include <fstream>
#include <iomanip>

#define _USE_MATH_DEFINES
#include <math.h>

void JumperObjective::init(const Eigen::VectorXd & x0)
{
  em = world->em_[0];
  conf->getFloat("beamwidth", opts.beamwidth);
  conf->getFloat("leftmargin", opts.leftmargin);
  materials.resize(em->m.size());
  for(unsigned int ii= 0; ii<em->m.size(); ii++){
    materials[ii] = (MaterialQuad*)em->m[ii];
  }
  state.fe = em->fe;
  state.fixedDof= em->fixedDof;
  state.X = em->X;
  state.lb=em->lb;
  if(em->dim== 2){
    segmentJumper(*em, opts);
  }else{
    segmentJumper3d(*em, opts);
  }
  setParam(x0);
  lambda.clear();
  lambda.resize(x0.size(), 0);
  g.resize(x0.size(), 0);
  rho.resize(x0.size(), 0);
  rho[0] = 1e3;
  rho[1] = 1e-3;
  updateG = true;
}

///@brief set simulation to initial. Restore external forces, x, v.
void JumperObjective::resetSim()
{
  em->fixedDof = state.fixedDof;
  em->lb = state.lb;
  em->fe = state.fe;
  em->x = em->X;
  std::fill(em->v.begin(), em->v.end(), Eigen::Vector3d::Zero());

  world->activeMesh = 0;
  world->frameCnt = 0;
  world->loading_t = 0;
  world->stage = World::SIM_LOADING;
  world->stepper = world->stepperList[0];
  world->stepperIdx = 0;
}

void JumperObjective::setParam(const Eigen::VectorXd & x0)
{
  param = x0;
  std::cout<<"JumperObjective::setParam \n";
  for(int ii=0; ii<x0.size();ii++){
    std::cout<<x0[ii]<<" ";
    opts.offset[ii] = x0[ii];
  }
  std::cout<<"\n";

  em->X = state.X;
  if(em->dim==2){
    stretchJumper(*em, opts);
  }else if(em->dim==3){
    stretchJumper3d(*em, opts);
  }
  em->initElements(materials[0]->q);
}

///@override
double JumperObjective::f()
{
  resetSim();
  
  std::string filename("opt_mesh");
  filename += std::to_string(nEval);
  filename += ".txt";
  std::ofstream out(filename);
  std::cout<<"save quad file "<<filename<<"\n";
  em->saveMesh(out);
  out.close();

  world->loop();

  double ang = world->totalAng;
  Eigen::Vector3d center = em->centerOfMass();
  bool addGravity = true;
  double Ek = em->getKineticEnergy(addGravity);
  
  //angular momentum weight
  //double wL = 2500;
  //double centerW = 1e2;
  //rotating clockwise.
  int NFlip = (int)(-ang / (2 * M_PI));
  double landAng = -ang - NFlip * 2 * M_PI;
  if (landAng > M_PI){
    landAng -= 2 * M_PI;
  }
  double val = Ek;
  double clearWallY = 0.1;
  double clearWallX = 0.01;
  std::vector<double> g0 = g;
  //hits wall
  //jumping onto
//  if(world->collideWith == 1 || std::abs(world->boxContactNormal[0]) > 0.1){
//    g[0] = clearWallY + world->box.mx[1] - center[1];
  //jumping over
  if (world->collideWith == 2 || center[0] <= world->box.mx[0] + clearWallX){
    g[0] = clearWallY + world->box.mx[1] - world->height;
  }
  else{
    g[0] = 0;
  }
  g[1] = 0;
  //g[1] = std::abs(landAng);
  
  //std::cout<<"final angle "<<ang <<" "<<landAng<<" "<<"\n";
  //std::cout << "Height " << world->height << "\n";
  //std::cout << "Energy " << Ek << "\n";
  //std::cout << "constriants " << g[0] << " " << g[1] << "\n";
  nEval ++;
  if (!updateG){
    g = g0;
  }
  val = landAng * landAng;
  for (int i = 0; i < g.size(); i++){
    double ghat = std::max(g[i] + lambda[i] / rho[i], 0.0);
    val += 0.5*rho[i] * ghat * ghat - 0.5 * lambda[i] * lambda[i] / rho[i];
  }

  return val;
}

Eigen::VectorXd JumperObjective::df()
{
  Eigen::VectorXd grad(param.size());
  Eigen::VectorXd x0 = param;
  updateG = false;
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
  updateG = true;
  return grad;
}
