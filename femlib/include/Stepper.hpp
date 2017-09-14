#ifndef STEPPER_HPP
#define STEPPER_HPP

#include <Eigen/Dense>

#include <thread>
#include <mutex>
#include <condition_variable>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

class ElementMesh;
struct Contact;
template<typename T> class LinSolver;

class Stepper{
public:

  Stepper();
  virtual void init();
  virtual int step();
  virtual void finalize(bool saveFile = false);
  virtual void resolveCollision(std::vector<Contact> & collision){}
  virtual int stepWrapper();

  ///@return positive if the time stepper converged.
  ///        negative if there is an error.
  ///        0 if needs to take more steps.
  virtual int oneStep() = 0;
  virtual ~Stepper();

  void logmsg(const std::string &msg);
  double totalAng();

  enum FrictionModel{INFINITE, NONE, TRESCA, FINITE};
  enum SimType{SIM_DYNAMIC, SIM_STATIC, SIM_RIGID};

  //=========================
  //simulation states
  ///@brief step length.
  float   h;

  ///@brief mesh for time stepping.
  ElementMesh * m;
  std::vector<ElementMesh * > meshes;
  ///@brief elapsed simulation time
  float   t0;
  float stiffness_t0;

  //=========================
  //simulation configurations

  bool handle_contact;

  ///@brief maximum number of time steps.
  int   nSteps;

  bool  addGravity;

  FrictionModel frictionModel;

  SimType  simType;

  LinSolver<double> *solver;

  std::vector < LinSolver<double> * > * solvers;
  ///@brief sequence of angles
  std::vector<float>           angSeq;

  ///@brief file name for a log file.
  std::string   logFileName;
  std::ofstream log;

};

#endif
