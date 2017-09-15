#include "StepperBDF2.hpp"
#include "ElementMesh.hpp"
#include "ArrayUtil.hpp"
#include "EigenUtil.hpp"
#include "femError.hpp"
#include "MaterialQuad.hpp"
#include "StrainEne.hpp"
#include "Timer.hpp"
#include "Eigen/Sparse"
#include "LinSolver.hpp"
#include "LinICPCG.hpp"
#include "LinPardiso.hpp"
#include "NewtonMin.hpp"
#include "Objective.hpp"
#include <sstream>

class BDF2Objective:public Objective
{
public:

  BDF2Objective():stepper(0){}

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
    grad = stepper->getGrad();
    if(fem_error<0){
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
    return stepper->m->dim * (int)stepper->m->X.size();
  }

  StepperBDF2 * stepper;
};

void StepperBDF2::init()
{
  Stepper::init();

  if (!m_Init){
    if(solver==0){
      std::cout<<"Warning. Using default cg linear solver\n";
      solver = new LinICPCG<double>();
    }
    solver->init();
    newton = new NewtonMin();
    BDF2Objective * obj = new BDF2Objective();
    newton->obj = obj;
    newton->solver = solver;
    newton->df_rel_tol = 1e-5f;
    newton->df_abs_tol = 1e-15f;
    newton->maxIters = 10;
    obj->stepper = this;
    m_Init = true;
  }

  D=Eigen::SparseMatrix<double>();
  p0=Eigen::VectorXd(m->dim * m->x.size());

  Eigen::SparseMatrix<double> K = m->getStiffnessSparse();
  solver->init(K);
  x0 = m->x;
  x1 = m->x;
  v0 = m->v;
  v1 = m->v;
}

StepperBDF2::~StepperBDF2()
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

StepperBDF2::StepperBDF2():newton(0)
{
  m_Init = false;
}

void StepperBDF2::updateStiffness()
{
  K = m->getStiffnessSparse();
  D = m->dampBeta * K + m->dampAlpha * m->M;
  K += (3/(2*h))*D;
  K = ((4.0/9.0)*h*h) * K + m->M;
}

int StepperBDF2::runNewton( std::vector<double> &x)
{
  std::stringstream msg;
  int newtonstatus = newton->minimize(x);
  int status = 0;
  if(newtonstatus<0){
    status = -1;
    msg.str("");
    msg<<"Newton collision err "<<newtonstatus<<"\n";
    logmsg(msg.str());
  }
  return status;
}

void StepperBDF2::update_p0()
{
  int dim = m->dim;
  //p0 = M(1/3 x1 - 1/3 x0 + h8/9v1-h2/9v0).
  Eigen::VectorXd vec(dim * m->X.size());
  for(int ii=0; ii<m->X.size(); ii++){
    for(int jj = 0; jj<dim; jj++){
      vec[ii*dim + jj] = (1.0/3.0)*x1[ii][jj] - (1.0/3.0) * x0[ii][jj] +
          h*(8.0/9.0)*v1[ii][jj] - h*(2.0/9.0)*v0[ii][jj];
    }
  }
  p0 = m->M*vec;
}

int StepperBDF2::oneStep()
{
  int dim = m->dim;
  int ndof = dim*(int)m->x.size();
  int status =0 ;

  updateStiffness();
  zeroOffDiag(K, m->fixedDof);
  x0 = x1;
  v0 = v1;
  x1 = m->x;
  v1 = m->v;
  update_p0();
  delta=Eigen::VectorXd::Zero(ndof);
  std::vector<double> x(ndof, 0);
  status = runNewton(x);

  //plot lagrange multipliers
  bool constrained = false;
  std::vector<double> grad = getGrad(constrained);
  for(unsigned int ii = 0; ii<m->x.size(); ii++){
    if(m->fixedDof[dim * ii+1]){
      for(int jj =0 ; jj<dim; jj++){
        double lambda = 2/h*grad[ii*dim + jj];
        m->fc[ii][jj] = lambda;
      }
    }else{
      m->fc[ii] = Eigen::Vector3d::Zero();
    }
  }

  return status;
}

float StepperBDF2::getEnergy()
{
  double E = 0;
  E  = (4.0/9.0)*h*h* m->getEnergy(addGravity);
  E -= delta.dot(p0);

  Eigen::VectorXd dx = Eigen::VectorXd::Map(delta.data(), delta.size());
  E += 0.5 * dx.dot(m->M * dx) + (1.0/3.0)*h*dx.dot(D * dx);

  return (float)E;
}

std::vector<double> StepperBDF2::getGrad(bool constrained)
{
  int dim = m->dim;
  int ndof = dim * (int)m->x.size();
  std::vector<double> rhs(ndof);
  std::vector<Eigen::Vector3d> F = m->getForce(addGravity);
  Eigen::VectorXd Dv, Mdx;
  Eigen::VectorXd vec(delta.size());
  Mdx = m->M * delta;

  for(unsigned int ii =0 ; ii<m->x.size(); ii++){
    for(int jj =0 ; jj<dim; jj++){
      vec[ii*dim + jj] = m->v[ii][jj];
    }
  }
  Dv = D*vec;
  //dV = Mdelta^i - p0 - 4/9h^2 F
  for(unsigned int ii =0 ; ii<F.size(); ii++){
    for(int jj =0 ; jj<dim; jj++){
      rhs[dim*ii+jj] = Mdx[dim*ii + jj]  - p0[ii*dim + jj]
          - (4.0/9.0)*h*h * (F[ii][jj] - Dv[ii*dim+jj]);
    }
  }
  if(constrained){
    setVal(rhs, m->fixedDof, 0.0);
  }
  return rhs;
}

Eigen::SparseMatrix<double>
StepperBDF2::getHessian()
{
  updateStiffness();
  Eigen::SparseMatrix<double> Kc = K;
  zeroOffDiag(Kc, m->fixedDof);
  return Kc;
}

void StepperBDF2::setx(double * x)
{
  int dim = m->dim;
  for(unsigned int ii =0; ii<m->x.size(); ii++){
    for(int jj =0; jj<dim; jj++){
      delta[ii*dim + jj] = x[ii*dim + jj];
      m->x[ii][jj] = x1[ii][jj] + delta[ii*dim + jj];
    }
    for(int jj =0 ; jj<dim; jj++){
      m->v[ii][jj] = (3.0/(2.0*h)) * (delta[ii*dim + jj] - (1.0/3.0)*(x1[ii][jj]-x0[ii][jj]));
    }
  }
}


void
StepperBDF2::resolveCollision(std::vector<bool> & collision)
{
  if(m->dim==2){
    resolveCollision2D(collision);
  }else if (m->dim==3){
    resolveCollision3D(collision);
  }
}

void
StepperBDF2::resolveCollision2D(std::vector<bool> & collision)
{
  
}

void
StepperBDF2::resolveCollision3D(std::vector<bool> & collision)
{
 
}
