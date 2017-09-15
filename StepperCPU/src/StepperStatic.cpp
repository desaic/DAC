#include "StepperStatic.hpp"
#include "Element.hpp"
#include "ElementMesh.hpp"
#include "ElementMeshUtil.hpp"
#include "Material.hpp"
#include "ArrayUtil.hpp"
#include "femError.hpp"
#include "LinSolver.hpp"
#include "LinPardiso.hpp"
#include "NewtonMin.hpp"
#include "Objective.hpp"
#include "EigenUtil.hpp"
#include "Objective.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <iostream>
#include <fstream>
#include <set>
#include <assert.h>

typedef Eigen::Triplet<double> Tripd;

class StaticObjective:public Objective
{
public:

  StaticObjective(int _n):n(_n),stepper(0){}

  void setx(double * x){
    stepper->setx(x);
  }

  double f(int & status)
  {
    status = 0;
    float val = stepper->getEnergy();
    if(fem_error != 0){
      fem_error = 0;
      status = -1;
    }
    return val;
  }

  int df(std::vector<double> & grad)
  {
    fem_error = 0;
    grad = stepper->getGrad();
    if(fem_error!=0){
      fem_error = 0;
      return -1;
    }
    return 0;
  }

  Eigen::SparseMatrix<double> d2f()
  {
    return stepper->getHessian();
  }

  int nx(){
    return n;
  }

  int n;
  StepperStatic * stepper;
};

StepperStatic::StepperStatic():
  force_norm0(-1),
  force_relative_tol(1e-5),
  force_abs_tol(1e-10)
{
  m_Init = false;
  simType = SIM_STATIC;
}

void StepperStatic::init()
{
  Stepper::init();
  int nConstraint = 0;
  for (int i = 0; i < CEq.size(); i++){
    nConstraint += CEq[i].size()-1;
  }
  int nVar = m->dim * (int)m->X.size() + nConstraint;

  if (!m_Init){
    if(solver==0){
      std::cout<<"StepperStatic Error: no linear solver is assigned.\n";
    }
    newton = new NewtonMin();

    StaticObjective * obj = new StaticObjective(nVar);
    newton->obj = obj;
    newton->solver = solver;
    newton->df_rel_tol = 1e-3f;
    newton->df_abs_tol = 1e-13f;
    newton->maxIters = 1;
    obj->stepper = this;
    m_Init = true;
  }

  force_norm0 = 1e-10f;
  std::vector<Eigen::Vector3d> F = m->getForce(addGravity);
  for (unsigned int ii = 0; ii<m->fe.size(); ii++){
    for (int jj = 0; jj<m->dim; jj++){
      if (m->fixedDof[m->dim * ii + jj]){
        continue;
      }
      force_norm0 = std::max(force_norm0, std::abs(F[ii][jj]));
    }
  }
}

StepperStatic::~StepperStatic()
{
  if(solver != 0){
    delete solver;
    solver = 0;
    delete newton->obj;
    newton->obj = 0;
    newton->solver = 0;
    delete newton;
  }
}

int StepperStatic::oneStep()
{
  int nConstraint = 0;
  for (int i = 0; i < CEq.size(); i++){
    nConstraint += CEq[i].size() - 1;
  }
  int nVar = m->dim * (int)m->X.size() + nConstraint;

  std::vector<double> x(nVar, 0);
  x0    = m->x;
  int newtonstatus = newton->minimize(x);
  //used to identify breaking vertices.
  std::vector<Eigen::Vector3d> F = m->getForce(addGravity);
  //std::cout << "static fixed0 size " << fixed0.size() << "\n";
  //std::cout << "=========handle contact " << handle_contact << "\n";
  int dim = m->dim;
  bool change=false;
  double forceInfNorm = 0;
  std::vector<Eigen::Vector3d> delta(x0.size());
  for(unsigned int ii=0; ii<F.size(); ii++){
    delta[ii] = m->x[ii] - x0[ii];
    m->fc[ii] = Eigen::Vector3d(0,0,0);
    for (int j = 0; j < m->dim; j++){
      if (m->fixedDof[dim*ii + j]>0){
        m->fc[ii] = F[ii];
      }
    }
    if(handle_contact && m->fixedDof[dim*ii+1] && F[ii][1]>0
      && (fixed0.size()<=dim*ii+1  || fixed0[dim*ii+1]==0)){
      for(int jj = 0; jj<dim; jj++){
        //std::cout << "remove contact " << ii << "\n";
        m->fixedDof[dim*ii+jj] = 0;
      }
    }else if( m->x[ii][1] < m->lb[ii][1] && m->fixedDof[dim*ii+1]==0
              && delta[ii][1]<0){
      change = true;
      m->x[ii][1] = m->lb[ii][1];
      for(int jj = 0; jj<dim; jj++){
        m->fixedDof[dim*ii+jj] = 1;
      }
    }
  }

  std::set<int> cset;
  for (int i = 0; i < CEq.size(); i++){
    for (int j = 0; j < CEq[i].size(); j++){
      cset.insert(CEq[i][j]);
    }
  }

  for(unsigned int ii=0; ii<F.size(); ii++){
    for(int jj = 0; jj<dim; jj++){
      auto c_it = cset.find(ii*dim + jj);
      if(! (m->fixedDof[ii*dim + jj] || c_it!=cset.end()) ){
        forceInfNorm = std::max(forceInfNorm, std::abs(F[ii][jj]));
      }
    }
  }
  double rel = forceInfNorm/force_norm0;
  std::cout << "rel " << rel << "\n";
  std::cout<<change<<" force inf norm "<<forceInfNorm<<" "<<force_norm0<<"\n";
  if( (!change)
    && (rel<force_relative_tol)){
    //if(!handle_contact){
    //  double eps = 1e-5;
    //  for(unsigned int ii = 0; ii<m->x.size(); ii++){
    //    if (m->x[ii][1] > m->lb[ii][1] + eps
    //      && (fixed0.size() <= dim*ii + 1 || fixed0[dim*ii + 1] == 0)){  //(!fixed0[dim*ii + 1]) && 
    //      for(int jj = 0; jj<dim; jj++){
    //        m->fixedDof[ii*dim + jj] = 0;
    //      }
    //    }
    //  }
    //}
    return 1;
  }
  return 0;
}

float
StepperStatic::getEnergy()
{
  float E =0 ;
  E = (float)m->getEnergy(addGravity);
  return E;
}

std::vector<double>
StepperStatic::getGrad(bool constrained)
{
  int dim = m->dim;
  int nConstraint = 0;
  for (int i = 0; i < CEq.size(); i++){
    nConstraint += CEq[i].size() - 1;
  }
  int nVar = dim * (int)m->X.size() + nConstraint;


  std::vector<double> rhs(nVar);
  std::vector<Eigen::Vector3d> F = m->getForce(addGravity);

  for(unsigned int ii =0 ; ii<F.size(); ii++){
    for(int jj =0 ; jj<dim; jj++){
      rhs[dim*ii+jj] = -F[ii][jj];
      if(constrained && m->fixedDof[dim*ii+jj]){
        rhs[ii*dim + jj] = 0;
      }
    }
  }
  Element * ele = m->e[0];
  double eleSize = m->x[ele->at(1)][dim-1] - m->x[ele->at(0)][dim-1];
  double tol = 0.01 * eleSize;
  if(handle_contact){
  for(unsigned int ii = 0; ii<m->x.size(); ii++){
    if(m->x[ii][1] < m->lb[ii][1] - tol){
      fem_error = -1;
//      std::cout<<ii<<"\n";
      break;
    }
  }
  }
  return rhs;
}

Eigen::SparseMatrix<double>
StepperStatic::getHessian()
{
  Eigen::SparseMatrix<double> K = m->getStiffnessSparse();
  zeroOffDiag(K, m->fixedDof);
  for (int i = 0; i < (int)CEq.size(); i++){
    equalityConstraint(CEq[i], m, K);
  }
  //double minDiag = 1e10;
  //for (int col=0; col<K.cols(); col++){
  //  //index of vertex at column col.
  //  for (Eigen::SparseMatrix<double>::InnerIterator it(K, col); it; ++it){
  //    int row = it.row();
  //    if(row==col){
  //      minDiag = std::min(it.valueRef(), minDiag);
  //    }
  //  }
  //}
  //std::cout << "Min Diagonal " << minDiag << "\n";

  return K;
}

void StepperStatic::setx(double *x)
{
  int dim = m->dim;

  for(unsigned int ii =0; ii<m->x.size(); ii++){
    for(int jj =0; jj<dim; jj++){
      m->x[ii][jj] = x0[ii][jj] + x[ii*dim + jj];
    }
  }
}
