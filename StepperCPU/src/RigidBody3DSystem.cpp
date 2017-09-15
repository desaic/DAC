// FlowableSystem.cpp
//
// Breannan Smith
// Last updated: 09/03/2015

#include "RigidBody3DSystem.h"
#include <iostream>
SparseMatrixsc formBodySpaceMassMatrix(const std::vector<RigidBody3DState*> & rigidBodies);
SparseMatrixsc formBodySpaceInverseMassMatrix(const std::vector<RigidBody3DState*> & rigidBodies);
SparseMatrixsc formWorldSpaceMassMatrix(const std::vector<RigidBody3DState*> & rigidBodies);
SparseMatrixsc formWorldSpaceInverseMassMatrix(const std::vector<RigidBody3DState*> & rigidBodies);

RigidBody3DSystem::~RigidBody3DSystem()
{}

///@brief assemble M0, Minv0
void RigidBody3DSystem::assemble()
{
  int nBodies = rigidBodies.size();
  int nqdof = 12 * nBodies;
  m_q.resize(nqdof);
  for(int ii = 0; ii<nBodies; ii++){
    m_q.segment<3>(3*ii) = rigidBodies[ii]->x;
    m_q.segment<9>(3*nBodies+9*ii) = Eigen::Map<VectorXs>(rigidBodies[ii]->R.data(), 9);
  }
  int nvdofs = 6 * nBodies;
  m_v.resize(nvdofs);
  for(int ii = 0; ii<nBodies; ii++){
    m_v.segment<3>(3*ii) = rigidBodies[ii]->v;
    m_v.segment<3>(3*nBodies + 3*ii) = rigidBodies[ii]->omega;
  }
  m_M0 = formBodySpaceMassMatrix(rigidBodies);
  m_Minv0 = formBodySpaceInverseMassMatrix(rigidBodies);
  m_M = formWorldSpaceMassMatrix(rigidBodies);
  m_Minv = formWorldSpaceInverseMassMatrix(rigidBodies);

  m_force.resize(nqdof);
}

///@brief disassemble states into struct.
void RigidBody3DSystem::setstate(const VectorXs & q, const VectorXs & v)
{
  int nBodies = rigidBodies.size();
  m_q = q;
  m_v = v;
  for(int body_idx = 0; body_idx<nBodies; body_idx++){
    rigidBodies[body_idx]->x = m_q.segment<3>(3*body_idx);
    rigidBodies[body_idx]->R = Eigen::Map<const Matrix33sr> (m_q.segment<9>(3*nBodies+9*body_idx).data());
  }
  for(int body_idx = 0; body_idx<nBodies; body_idx++){
    rigidBodies[body_idx]->v = m_v.segment<3>(3*body_idx);
    rigidBodies[body_idx]->omega = m_v.segment<3>(3*nBodies + 3*body_idx);
  }

  //only need to update world space mass;
  m_M = formWorldSpaceMassMatrix(rigidBodies);
  m_Minv = formWorldSpaceInverseMassMatrix(rigidBodies);
}

// Given positions and velocities, computes the force acting on the DoFs, overwriting the contents of F
void RigidBody3DSystem::computeForce( const VectorXs& q, const VectorXs& v, const scalar& t, VectorXs& F )
{
  int nBodies = rigidBodies.size();
  m_force.resize(v.size());
  m_force.setZero(m_force.size());
  for(int body_idx = 0; body_idx < nBodies; body_idx++){
    scalar M = rigidBodies[body_idx]->M;
    m_force.segment<3>(3*body_idx) = M*m_g + rigidBodies[body_idx]->f_linear;
  }
  F = m_force;
}

// Returns the mass matrix as a sparse matrix
const SparseMatrixsc& RigidBody3DSystem::M() const
{
  return m_M;
}

// Returns the inverse mass matrix as a sparse matrix
const SparseMatrixsc& RigidBody3DSystem::Minv() const
{
  return m_Minv;
}

// Returns the reference configuration mass matrix as a sparse matrix
const SparseMatrixsc& RigidBody3DSystem::M0() const
{
  return m_M0;
}

// Returns the reference configuration inverse mass matrix as a sparse matrix
const SparseMatrixsc& RigidBody3DSystem::Minv0() const
{
  return m_Minv0;
}

//static functions -------------------
SparseMatrixsc formBodySpaceMassMatrix(const std::vector<RigidBody3DState *> & rigidBodies)
{
  //body space mass
  int nvdofs = 6*rigidBodies.size();
  int nBodies = rigidBodies.size();
  SparseMatrixsc M0(nvdofs, nvdofs);
  M0.reserve(VectorXi::Ones(nvdofs));
  for(int body_idx = 0; body_idx<nBodies; body_idx ++){
    scalar M = rigidBodies[body_idx]->M;
    for(int dof_idx = 0; dof_idx<3; dof_idx++){
      int col = 3*body_idx + dof_idx;
      int row = col;
      M0.insert(row, col)=M;
    }
  }

  for(int body_idx = 0; body_idx<nBodies; body_idx ++){
    Vector3s I0 = rigidBodies[body_idx]->I0;
    for(int dof_idx = 0; dof_idx<3; dof_idx++){
      int col = 3*nBodies + 3*body_idx + dof_idx;
      int row = col;
      M0.insert(row, col)=I0[dof_idx];
    }
  }
  M0.makeCompressed();
  return M0;
}

SparseMatrixsc formBodySpaceInverseMassMatrix(const std::vector<RigidBody3DState *> &rigidBodies)
{
  //body space mass
  int nvdofs = 6*rigidBodies.size();
  int nBodies = rigidBodies.size();
  SparseMatrixsc M0(nvdofs, nvdofs);
  M0.reserve(VectorXi::Ones(nvdofs));
  for(int body_idx = 0; body_idx<nBodies; body_idx ++){
    scalar M = 1.0/rigidBodies[body_idx]->M;
    for(int dof_idx = 0; dof_idx<3; dof_idx++){
      int col = 3*body_idx + dof_idx;
      int row = col;
      M0.insert(row, col)=M;
    }
  }

  for(int body_idx = 0; body_idx<nBodies; body_idx ++){
    Vector3s I0 = rigidBodies[body_idx]->I0;
    for(int dof_idx = 0; dof_idx<3; dof_idx++){
      int col = 3*nBodies + 3*body_idx + dof_idx;
      int row = col;
      M0.insert(row, col) = 1.0/I0[dof_idx];
    }
  }
  M0.makeCompressed();
  return M0;
}

SparseMatrixsc formWorldSpaceMassMatrix(const std::vector<RigidBody3DState*> & rigidBodies)
{
  int nvdofs = 6*rigidBodies.size();
  int nBodies = rigidBodies.size();
  SparseMatrixsc Mbody(nvdofs, nvdofs);
  {
    VectorXi column_nonzeros{ nvdofs };
    column_nonzeros.segment( 0, 3 * nBodies ).setOnes();
    column_nonzeros.segment( 3 * nBodies, 3 * nBodies ).setConstant( 3 );
    Mbody.reserve( column_nonzeros );
  }
  for(int body_idx = 0; body_idx<nBodies; body_idx ++){
    scalar M = rigidBodies[body_idx]->M;
    for(int dof_idx = 0; dof_idx<3; dof_idx++){
      int col = 3*body_idx + dof_idx;
      int row = col;
      Mbody.insert(row, col)=M;
    }
  }

  for(int body_idx = 0; body_idx<nBodies; body_idx ++){
    Matrix33sr Rmat = rigidBodies[body_idx]->R;
    Vector3s I0 = rigidBodies[body_idx]->I0;
    Matrix33sr I = Rmat * I0.asDiagonal() * Rmat.transpose();
    for( unsigned row_idx = 0; row_idx < 3; ++row_idx ){
      const unsigned col{ 3 * nBodies + 3 * body_idx + row_idx };
      for( unsigned col_idx = 0; col_idx < 3; ++col_idx ){
        const unsigned row{ 3 * nBodies + 3 * body_idx + col_idx };
        Mbody.insert(row, col) = I(row_idx, col_idx);
      }
    }
  }
  Mbody.makeCompressed();
  return Mbody;
}

SparseMatrixsc formWorldSpaceInverseMassMatrix(const std::vector<RigidBody3DState*> & rigidBodies)
{
  int nvdofs = 6*rigidBodies.size();
  int nBodies = rigidBodies.size();
  SparseMatrixsc Mbody(nvdofs, nvdofs);
  {
    VectorXi column_nonzeros{ nvdofs };
    column_nonzeros.segment( 0, 3 * nBodies ).setOnes();
    column_nonzeros.segment( 3 * nBodies, 3 * nBodies ).setConstant( 3 );
    Mbody.reserve( column_nonzeros );
  }
  for(int body_idx = 0; body_idx<nBodies; body_idx ++){
    scalar M = 1.0/rigidBodies[body_idx]->M;
    for(int dof_idx = 0; dof_idx<3; dof_idx++){
      int col = 3*body_idx + dof_idx;
      int row = col;
      Mbody.insert(row, col)=M;
    }
  }

  for(int body_idx = 0; body_idx<nBodies; body_idx ++){
    Matrix33sr Rmat = rigidBodies[body_idx]->R;
    Vector3s I0 = rigidBodies[body_idx]->I0;
    Matrix33sr Iinv = Rmat * I0.array().inverse().matrix().asDiagonal() * Rmat.transpose();
    for( unsigned row_idx = 0; row_idx < 3; ++row_idx ){
      const unsigned col{ 3 * nBodies + 3 * body_idx + row_idx };
      for( unsigned col_idx = 0; col_idx < 3; ++col_idx ){
        const unsigned row{ 3 * nBodies + 3 * body_idx + col_idx };
        Mbody.insert(row, col) = Iinv(row_idx, col_idx);
      }
    }
  }
  Mbody.makeCompressed();
  return Mbody;
}

void RigidBody3DSystem::printState(){
  std::cout<<"m_q"<<m_q<<"\n";
  std::cout<<"m_v"<<m_v<<"\n";
  std::cout<<"M0"<<m_M0<<"\n";
  std::cout<<"Minv0"<<m_Minv0<<"\n";
  std::cout<<"M"<<m_M<<"\n";
  std::cout<<"Minv"<<m_Minv<<"\n";
}
