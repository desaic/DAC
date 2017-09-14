#include "Stepper.hpp"
#include "Timer.hpp"
#include "ElementMesh.hpp"
#include "ArrayUtil.hpp"
#include "EigenUtil.hpp"
#include "FileUtil.hpp"

#include <atomic>
#include <iostream>
#include <fstream>

# define M_PI           3.14159265358979323846

Stepper::Stepper() :m(0),
  handle_contact(false),
  nSteps(100), h(1.0f),
  t0(0), stiffness_t0(0),
  addGravity(false),
  frictionModel(INFINITE),
  simType(SIM_DYNAMIC),
  solvers(0)
{}

Stepper::~Stepper(){}

void Stepper::init()
{
  t0=0;
  stiffness_t0 = 0;
  angSeq.clear();
  log.open(logFileName);
//  if(!log.good()){
//    std::cout<<"cannot open log file "<<logFileName<<"\n";
//  }
}

int Stepper::stepWrapper()
{
  if(m==0){
    std::cout<<"Need a mesh in stepper.\n";
    return -1;
  }
  int ret = 0;

  Timer t;
  t.start();
  if(angSeq.size() == 0 && m->topVerts.size()>0){
    double ang = m->getTopAng();
    angSeq.push_back((float)ang);
  }
  ret = oneStep();
  t.end();
  float duration = t.getSeconds();
//  std::cout<<" sim time: "<<t0<<
//             " computation time: " << duration << "\n";

  if(m->topVerts.size()>0){
    double ang = m->getTopAng();
    angSeq.push_back((float)ang);
//    std::cout<<"angle : " << ang<<"=================\n";
  }
  return ret;
}

int Stepper::step()
{
  int ret = 0;

  for (int ii = 0; ii < nSteps; ii++){
    ret = stepWrapper();
    if (ret < 0){
//      break;
    }
    if (ret>0){
  //    std::cout << "Early Termination." << std::endl;
      break;
    }
//    std::cout<<m->getTopAng()<<"\n";
  }

//  std::cout << "End stepper loop\n";
//  finalize();
//  std::cout << "Finish stepper\n";
  return ret;
}


void Stepper::logmsg(const std::string & msg)
{
  if(log.is_open()){
    log<<msg<<"\n";
  }else{
    std::cout<<"cannot write to log\n";
  }
}

double Stepper::totalAng()
{
  float maxChange = 3;
  double totalAng = 0;
  if(angSeq.size()>0){
//    std::cout<<angSeq[0]<<"\n";
    for(unsigned int ii = 1; ii<angSeq.size(); ii++){
//      std::cout<<angSeq[ii]<<"\n";
      double dtheta = angSeq[ii] - angSeq[ii-1];
      if(dtheta>maxChange){
        dtheta -= 2*M_PI;
      }
      if(dtheta<-maxChange){
        dtheta += 2*M_PI;
      }
      totalAng += dtheta;
    }
  }
  return totalAng;
}

void Stepper::finalize(bool saveFile)
{}
