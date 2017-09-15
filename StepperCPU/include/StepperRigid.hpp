#ifndef STEPPERRIGID_HPP
#define STEPPERRIGID_HPP

#include "Stepper.hpp"
#include "DMVMap.h"
#include "RigidBody3DSystem.h"

#include <vector>

class StepperRigid:public Stepper{
public:
  StepperRigid();
  ~StepperRigid();
  void init();
  int oneStep();

  DMVMap dmv;
  RigidBody3DSystem fsys;
  RigidBody3DState rigid;
};

#endif
