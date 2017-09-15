#include "StepperLinNewmark.hpp"

#include "ElementMesh.hpp"
#include "ArrayUtil.hpp"
#include "EigenUtil.hpp"
#include "femError.hpp"

#include "Timer.hpp"
#include "Eigen/Sparse"
#include "LinSolver.hpp"
#include "LinICPCG.hpp"
#include "LinPardiso.hpp"
#include <fstream>
#include <assert.h>

StepperLinNewmark::StepperLinNewmark()
{
  h=0.001f;
  m_Init = false;
}

int StepperLinNewmark::oneStep()
{
  int dim = m->dim;
  int ndof = dim*(int)m->x.size();

  //acceleration at step n.
  std::vector<double> a_n(ndof);

  std::vector<Eigen::Vector3d> force;

  force = m->getForce(addGravity);

  if(anext.size() == a_n.size()){
    a_n = anext;
  }else{
    solve_acceleration(force, std::vector<bool>(0), a_n);
  }

  for(unsigned int ii =0; ii<m->x.size(); ii++){
    for(int jj =0 ; jj<dim; jj++){
      m->x[ii][jj] += 0.5*h*h*a_n[dim*ii+jj];
    }
    m->x[ii] += h*m->v[ii];
  }
  anext.resize(a_n.size(), 0);
  force = m->getForce(addGravity);
  solve_acceleration(force, std::vector<bool>(0), anext);

  for(unsigned int ii = 0; ii<force.size(); ii++){
    for(int jj =0 ; jj<dim; jj++){
      m->v[ii][jj] += 0.5*h*(a_n[dim*ii+jj] + anext[dim*ii+jj]);
    }
  }

  return 0;
}

double StepperLinNewmark::solve_acceleration(const std::vector<Eigen::Vector3d> & forces,
                                          const std::vector<bool> & collide,
                                          std::vector<double> &bb)
{
  if (!m_Init){
    if(solver==0){
      std::cout<<"Warning. Using default cg linear solver\n";
      solver = new LinICPCG<double>();
    }
    solver->init();
    m_Init = true;
  }

  float beta = 0.25;
  int dim = m->dim;
  int ndof = dim * (int)m->x.size();
  std::vector<double> rhs(ndof);
  std::vector<double> x(ndof,0);

  assert(m && ndof==bb.size());

  Eigen::SparseMatrix<double> K = (m->getStiffnessSparse()).cast<double>();

  for(unsigned int ii = 0;ii<m->x.size(); ii++){
    for(int jj = 0;jj<dim;jj++){
      int row = dim*ii + jj;
      rhs[ row ] = forces[ii][jj];
    }
  }
  setVal(rhs, m->fixedDof, 0.0);
  K *= beta * h * h;
  K += m->M;
  zeroOffDiag(K, m->fixedDof);
  double residual;
  int iters;
  int ret = solver->solve(K, rhs, x, residual, iters);

  for(int ii = 0;ii<x.size();ii++){
    bb[ii] = x[ii];
  }

  return 0;
} 
