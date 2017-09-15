#include "StepperRigid.hpp" 
#include "ElementMesh.hpp"
#include "EigenUtil.hpp"
#include "MomentTools.h"
#include "TrigMesh.hpp"

void computeMoments(ElementMesh *m, Vector3s & I, Vector3s& center, Matrix33sr& R );

void getFrame(const ElementMesh & m, RigidBody3DState & rigid);

StepperRigid::StepperRigid()
{
  simType = SIM_RIGID;
}

StepperRigid::~StepperRigid()
{}

void StepperRigid::init()
{
//  if(m->dim == 2){
//    std::cout<<"Error in StepperRigid. Does not work for 2D\n";
//    return;
//  }
  Stepper::init();
  std::vector<Eigen::Vector3d> x = m->x;
  m->x = m->X;
  rigid = RigidBody3DState();
  fsys = RigidBody3DSystem();
  computeMoments(m, rigid.I0, rigid.x0, rigid.R0);
  for(unsigned int ii = 0; ii<m->v.size(); ii++){
    m->X[ii] -= rigid.x0;
  }
  m->x = x;
  getFrame(*m, rigid);
  fsys.addRigidBody(&rigid);
  fsys.assemble();
  rigid.transform(m->X, m->x);
}

int StepperRigid::oneStep()
{
  VectorXs q0, v0,q1, v1;
  q0 = fsys.q();
  v0 = fsys.v();
  dmv.flow(q0, v0, fsys, 1, h, q1, v1);
  fsys.setstate(q1,v1);
  rigid.transform(m->X, m->x);
  for(unsigned int ii = 0; ii<m->v.size(); ii++){
    for(int jj = 0; jj<m->dim;jj++){
      m->v[ii][jj] = v1[jj];
    }
  }
  Eigen::Vector3d center = m->centerOfMass();
  for(unsigned int ii = 0; ii<m->v.size(); ii++){
    Eigen::Vector3d r = m->x[ii] - center;
    m->v[ii] += v1.segment<3>(3).cross(r);
  }

  return 0;
}

void getFrame(const ElementMesh & m, RigidBody3DState & rigid)
{
  Eigen::Matrix3Xd X, x;
  X=assemble3n(m.X);
  x=assemble3n(m.x);
  Eigen::Matrix3d R;
  if(m.dim == 2){
    R = shapeMatch2d(x,X);
  }else{
    R = shapeMatch(x,X);
  }
  Eigen::Vector3d center = m.centerOfMass();
  Eigen::Vector3d p = m.linearMomentum();

  rigid.v = p/m.totalMass;
  //std::cout << " v " << rigid.v << "\n";
  rigid.M = m.totalMass;
  rigid.x = center;
  rigid.R = R.transpose() * rigid.R0;

  Matrix33sr Iinv = rigid.R * rigid.I0.array().inverse().matrix().asDiagonal() * rigid.R.transpose();
  Eigen::Vector3d L = m.angularMomentum(center);
  rigid.omega= 0.98*Iinv * L;
}

void computeMoments(ElementMesh * m, Vector3s & I, Vector3s& center, Matrix33sr &R )
{
  Matrix3s Rarg;
  MomentTools::computeMoments( m, I, center, Rarg);
  R = Rarg;
}
