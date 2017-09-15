#include "StepperNewtonDyn.hpp"
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

StepperNewtonDyn::StepperNewtonDyn():dx_tol(1e-5f)
{
  h=0.001f;
  m_Init = false;
}

int StepperNewtonDyn::oneStep()
{
    std::vector<Eigen::Vector3d> force = m->getForce(addGravity);

  int dim = m->dim;
  int ndof = dim*(int)m->x.size();
  std::vector<double> bb(ndof);

  for(unsigned int ii =0 ; ii<force.size(); ii++){
    force[ii] = h*h*force[ii] + h*m->mass[ii]*m->v[ii];
  }

  compute_dx_sparse(m, force, std::vector<bool>(0), bb);
  std::vector<bool> collide(m->x.size(),false);
  bool hasCollision = false;
  //hard-coded collision detection.
  for(unsigned int ii = 0; ii<m->x.size(); ii++){
    if( bb[dim*ii+1] + m->x[ii][1] < 0 ){
//      std::cout<<ii<<"\n";
      collide[ii] = true;
      hasCollision = true;
    }
  }
  if(hasCollision){
    compute_dx_sparse(m, force, collide, bb);
    //force sanity check
//    force = m->getForce();
//    for(unsigned int ii =0; ii<m->x.size(); ii++){
//      if(collide[ii] && force[ii][1]>0){
//        std::cout<<"Force check "<<ii<<"\n";
//      }
//    }
  }

  for(unsigned int ii =0; ii<force.size(); ii++){
    for(int jj =0 ; jj<dim; jj++){
      m->x[ii][jj] += bb[dim*ii+jj];
      m->v[ii][jj] = (1.0/h) * bb[dim*ii+jj];
    }
  }
  return 0;
}

float StepperNewtonDyn::compute_dx_sparse(ElementMesh * iMesh,
                                          const std::vector<Eigen::Vector3d> &iForces,
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

  int dim = iMesh->dim;
  int ndof = dim * (int)iMesh->x.size();

  assert(iMesh && ndof==bb.size());

  Timer timer;
  timer.start();
  Eigen::SparseMatrix<double> K = (iMesh->getStiffnessSparse()).cast<double>();
  timer.end();
//  std::cout<<"Assembly time "<< timer.getSeconds()<<"\n";

  timer.start();

  for(unsigned int ii = 0;ii<m->x.size(); ii++){
    for(int jj = 0;jj<dim;jj++){
      int row = dim*ii + jj;
      bb[ row ] = iForces[ii][jj];
    }
  }

  setVal(bb, m->fixedDof, 0.0);
  int nrow = ndof;
  int ncol = ndof;

  std::vector<double> rhs(nrow);
  std::vector<double> x(nrow,0);
  K *= h*h;
  K = K + m->M;
  zeroOffDiag(K, m->fixedDof);
  timer.end();
//  std::cout<<"Copy data time "<< timer.getSeconds()<<"\n";
  timer.start();

  //collision correction
  if(collide.size() == iMesh->x.size()){
    for (int col=0; col<ncol;col++){
      int vcol = col/dim;
      if(collide[vcol]){
        rhs[col] = 0;
      }
      for (Eigen::SparseMatrix<double>::InnerIterator it(K, col); it; ++it){
       int row = it.row();
       int vrow = row/dim;
       if(collide[vrow] || collide [vcol]){
         if(row != col){
           it.valueRef()=0;
         }else{
//           std::cout<<mat.value[col][ i-m_I[col] ]<<"\n";
         }
       }
      }
    }
  }

  timer.end();
//  std::cout<<"Copy collision time "<< timer.getSeconds()<<"\n";

  double residual=0;
  int iters=0;
  timer.start();
  int ret = solver->solve(K, rhs, x, residual, iters);
//  std::cout<< "lin solver ret " <<ret<<" iters " <<
//              iters << " residual "<<residual<<"\n";
  timer.end();
//  std::cout<<"Solver time "<< timer.getSeconds()<<"\n";

  for(int ii = 0;ii<x.size();ii++){
    bb[ii] = x[ii];
  }

  return 0;
}
