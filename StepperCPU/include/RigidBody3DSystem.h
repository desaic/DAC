// FlowableSystem.h
//
// Breannan Smith
// Last updated: 09/03/2015

#ifndef RIGIDBODY3D_SYSTEM_H
#define RIGIDBODY3D_SYSTEM_H

#include "MathDefines.h"

struct RigidBody3DState
{

  RigidBody3DState():x(0,0,0),v(0,0,0),omega(0,0,0),f_linear(0,0,0),tau(0,0,0),
  R(Matrix3s::Identity())
  {}
  //quantities do not evolve with simulation

  ///@brief body space diagonal inertia
  Vector3s I0;

  ///@brief rotation from body space to rest position of the mesh at t0.
  Matrix33sr R0;

  ///@brief mass.
  scalar M;

  Vector3s x0;

  //quantities that change with simulation
  ///@brief center of mass
  Vector3s x;
  ///@brief rotation from body space to current simulation state. Row major.
  Matrix33sr R;

  Vector3s v;

  Vector3s omega;

  ///@brief force at center of mass
  Vector3s f_linear;

  ///@brief torque
  Vector3s tau;

  ///@param vertices0 assume in body space centered at 0
  void transform(const std::vector<Vector3s> & vertices0, std::vector<Vector3s> & vertices)
  {
    assert(vertices0.size() == vertices.size());
    Matrix3s Rprod = R * R0.transpose();
    for(unsigned int ii=0; ii<vertices.size(); ii++){
      vertices[ii] = Rprod*vertices0[ii] + x;
    }
  }
};

class RigidBody3DSystem
{

public:

  RigidBody3DSystem():m_g(0,-9.8,0){}
  virtual ~RigidBody3DSystem();

  // Given positions and velocities, computes the force acting on the DoFs, overwriting the contents of F
  void computeForce( const VectorXs& q, const VectorXs& v, const scalar& t, VectorXs& F );

  ///@brief assemble M0, Minv0, q, v
  void assemble();

  ///@brief Disassemble global position and velocity. Update M Minv.
  void setstate(const VectorXs &q, const VectorXs &v);

  // Returns the mass matrix as a sparse matrix
  const SparseMatrixsc& M() const ;
  // Returns the inverse mass matrix as a sparse matrix
  const SparseMatrixsc& Minv() const ;

  // Returns the reference configuration mass matrix as a sparse matrix
  const SparseMatrixsc& M0() const ;
  // Returns the reference configuration inverse mass matrix as a sparse matrix
  const SparseMatrixsc& Minv0() const ;

  const VectorXs q()const {
    return m_q;
  }

  const VectorXs v()const{
    return m_v;
  }

  void addRigidBody(RigidBody3DState * r){
    rigidBodies.push_back(r);
  }

  RigidBody3DState * getRigidBody(int idx){
    assert(idx<rigidBodies.size());
    return rigidBodies[idx];
  }

  void printState();
  Vector3s m_g;

private:

  VectorXs m_force;
  ///@brief positions of all bodies followed by rotation matrices
  VectorXs m_q;
  VectorXs m_v;

  SparseMatrixsc m_M0;
  SparseMatrixsc m_Minv0;
  SparseMatrixsc m_M;
  SparseMatrixsc m_Minv;

  std::vector<RigidBody3DState*> rigidBodies;
};

#endif
