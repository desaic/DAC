#include "StepperNewmark.hpp"
#include "ElementMesh.hpp"
#include "ArrayUtil.hpp"
#include "EigenUtil.hpp"
#include "Contact.hpp"
#include "femError.hpp"
#include "MaterialQuad.hpp"
#include "StrainEne.hpp"
#include "Timer.hpp"
#include "Eigen/Sparse"
#include "LinSolver.hpp"
#include "NewtonMin.hpp"
#include "Objective.hpp"
#include <sstream>
#include <cstdint>
class NewmarkObjective:public Objective
{
public:

  NewmarkObjective():stepper(0){}

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

  StepperNewmark * stepper;
};

void StepperNewmark::init()
{
  Stepper::init();
  if (!m_Init){
    newton = new NewtonMin();
    NewmarkObjective * obj = new NewmarkObjective();
    newton->obj = obj;
    newton->df_rel_tol = 1e-5f;
    newton->df_abs_tol = 1e-15f;
    newton->maxIters = 10;
    obj->stepper = this;
    m_Init = true;
  }
}

StepperNewmark::~StepperNewmark()
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

StepperNewmark::StepperNewmark() :newton(0), useOldContact(false)
{
  m_Init = false;
}

void StepperNewmark::updateStiffness()
{
  std::vector <float> dcoef(m->m.size(), 1e-3);
  dcoef[0] = m->dampBeta;
  if (dcoef.size() > 1){
    dcoef[1] = 20 * m->dampBeta;
  }
  //K = m->getStiffnessSparse();
  //DampMat = m->dampBeta * K + m->dampAlpha * m->M;
  K = m->getStiffnessDamping(DampMat, dcoef);
//  Timer timer;
//  timer.startWall();
  K += (2 / h)*DampMat;
  K = (0.25*h*h) * K + m->M;
//  timer.endWall();
//  std::cout<<"newmark Eigen assemble time "<<timer.getSecondsWall()<<"\n";
}

int StepperNewmark::runNewton( std::vector<double> &x)
{
  std::stringstream msg;
  int newtonstatus = newton->minimize(x);
  return newtonstatus;
}

void StepperNewmark::updateMb()
{
  std::vector<Eigen::Vector3d> force;
  int dim = m->dim;
  int ndof = dim*(int)m->x.size();
  Eigen::VectorXd v(ndof), Dv(ndof);
  force = m->getForce(addGravity);

  for(unsigned int ii = 0; ii<m->v.size(); ii++){
    for(int jj =0; jj<dim; jj++){
      v[ii*dim + jj] = m->v[ii][jj];
    }
  }

  Dv = h * (m->M * v);
  if (m->vu.size() > 0){
    std::cout << "use unprojected v for damping.\n";
    for (unsigned int ii = 0; ii<m->v.size(); ii++){
      for (int jj = 0; jj<dim; jj++){
        v[ii*dim + jj] = m->vu[ii][jj];
      }
    }
  }
  Dv -= (0.25*h*h) * (DampMat*v);
  Mb.resize(force.size());
  for(unsigned ii=0; ii<Mb.size(); ii++){
    Mb[ii] = 0.25*h*h*force[ii];
    for(int jj =0 ; jj<dim; jj++){
      Mb[ii][jj] += Dv[ii*dim + jj];
    }
  }
}

void checkDamping(ElementMesh * em, StepperNewmark * stepper)
{
  em->x[0][0] += 0.01;
  stepper->updateStiffness();
  for(unsigned int ii = 0; ii<em->v.size(); ii++){
    em->v[ii] = Eigen::Vector3d(ii,0.5*ii,0);
  }
  int dim = em->dim;
  Eigen::VectorXd v(dim * em->v.size());
  for(unsigned int ii = 0; ii<em->v.size(); ii++){
    v.block(dim*ii,0,dim,1) = em->v[ii].block(0,0,dim,1);
  }
  Eigen::VectorXd f = stepper->DampMat*v;
  Eigen::MatrixXd R = Eigen::MatrixXd(3,3);
  double theta = 0.1;
  R<<cos(theta), -sin(theta),0,sin(theta), cos(theta),0,0,0,1;
  for(unsigned int ii=0; ii<em->x.size(); ii++){
    em->x[ii] = R*em->x[ii];
  }
  stepper->updateStiffness();
  for(unsigned int ii = 0; ii<em->v.size(); ii++){
    Eigen::Vector3d vi = (R*em->v[ii]);
    v.block(dim*ii,0,dim,1) = vi.block(0,0,dim,1);
  }
  Eigen::VectorXd f1 = stepper->DampMat*v;
  for(unsigned int ii = 0; ii<em->v.size(); ii++){
    Eigen::Vector3d vi(0,0,0);
    vi.block(0,0,dim,1) = f1.block(dim*ii,0,dim,1);
    vi = (R.transpose()*vi);
    f1.block(dim*ii,0,dim,1) = vi.block(0,0,dim,1);
  }
  for(int ii= 0 ; ii<f1.rows(); ii++){
    std::cout<<f[ii]<<" "<<f1[ii]<<"\n";
  }
}

int StepperNewmark::oneStep()
{
  int dim = m->dim;
  int ndof = dim*(int)m->x.size();
  int status =0 ;

  if (solver == 0){
    std::cout << "Newmark need a linear solver\n";
    return -1;
  }
  newton->solver = solver;
  updateStiffness();
  zeroOffDiag(K, m->fixedDof);
  updateMb();
  x0 = m->x;
  v0 = m->v;

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

float StepperNewmark::getEnergy()
{
  //not implemented.
  double E =0;
  return (float)E;
}

std::vector<double> StepperNewmark::getGrad(bool constrained)
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
  Dv = DampMat*vec;
  //dV = Mq^i - Mb - 1/4h^2 F
  for(unsigned int ii =0 ; ii<F.size(); ii++){
    for(int jj =0 ; jj<dim; jj++){
      rhs[dim*ii+jj] = Mdx[dim*ii + jj]  - Mb[ii][jj]
          - 0.25*h*h * (F[ii][jj] - Dv[ii*dim+jj]);
    }
  }
  if(constrained){
    setVal(rhs, m->fixedDof, 0.0);
  }
  return rhs;
}

std::vector<double> StepperNewmark::getRHS()
{
  int dim = m->dim;
  int ndof = dim * (int)m->x.size();
  std::vector<double> rhs(ndof);
  std::vector<Eigen::Vector3d> F = m->getForce(addGravity);
  Eigen::VectorXd Dv;
  Eigen::VectorXd vec(delta.size());

  for (unsigned int ii = 0; ii<m->x.size(); ii++){
    for (int jj = 0; jj<dim; jj++){
      vec[ii*dim + jj] = m->v[ii][jj];
    }
  }
  Dv = DampMat*vec;
  //dV = Mq^i - Mb - 1/4h^2 F
  for (unsigned int ii = 0; ii<F.size(); ii++){
    for (int jj = 0; jj<dim; jj++){
      rhs[dim*ii + jj] = Mb[ii][jj] + 0.25*h*h * (F[ii][jj] - Dv[ii*dim + jj]);
    }
  }
  setVal(rhs, m->fixedDof, 0.0);
  return rhs;
}


Eigen::SparseMatrix<double>
StepperNewmark::getHessian()
{
  updateStiffness();
  Eigen::SparseMatrix<double> Kc = K;
  zeroOffDiag(Kc, m->fixedDof);
  return Kc;
}

void StepperNewmark::setx(double * x)
{
  int dim = m->dim;
  for(unsigned int ii =0; ii<m->x.size(); ii++){
    for(int jj =0; jj<dim; jj++){
      delta[ii*dim + jj] = x[ii*dim + jj];
      m->x[ii][jj] = x0[ii][jj] + delta[ii*dim + jj];
    }
    for(int jj =0 ; jj<dim; jj++){
      m->v[ii][jj] = (2.0/h) * delta[ii*dim + jj] - v0[ii][jj];
    }
  }
}


void
StepperNewmark::resolveCollision(std::vector<Contact> & collision)
{
  //if(m->dim==2){
  //  resolveCollision2D(collision);
  //}else if (m->dim==3){
    resolveCollision3D(collision);
  //}
}

void ortho3(const Eigen::Vector3d & x,
  Eigen::Vector3d & y,
  Eigen::Vector3d & z)
{
  y = Eigen::Vector3d(0, 1, 0);
  if (std::abs(x[1]) > 0.9){
    y = Eigen::Vector3d(1, 0, 0);
  }
  z = x.cross(y);
#ifdef _DEBUG //check z has a length
  if(z.squaredNorm()<1e-20){
    std::cout << "Computing frame: axis has short length.\n";
  }
#endif 
  z.normalize();
  y = z.cross(x);
  y.normalize();
}

void printVec(double * a, int len)
{
	for (int i = 0; i < len; i++) {
		std::cout << a[i] << " ";
		if (i % 3 == 2) {
			std::cout << "\n";
		}
	}
}

int
StepperNewmark::resolveCollision3D_old(std::vector<Contact> & collision)
{
  int dim = m->dim;
  int N_NEWTON = 4;
  std::vector<Eigen::SparseMatrix<double> >DampMats(meshes.size()), Ks(meshes.size());
  std::vector<std::vector<Eigen::Vector3d> > Mbs(meshes.size());
  //std::cout << collision.size() << " contacts.\n";
  //for (size_t i = 0; i < collision.size(); i++){
  //  std::cout << collision[i].v1 << "\n";
  //}
  //number of constraints
  int NC = (int)collision.size();
  //total dof
  int nDof = 0;
  for (int mi = 0; mi < meshes.size(); mi++){
    m = meshes[mi];
    nDof += dim * (int)m->X.size();
  }
  //back solve to compute Ntilde and Dtilde
  Eigen::MatrixXd Ntilde(nDof, NC), Dtilde[2];
  Dtilde[0] = Eigen::MatrixXd(nDof, NC); //concatenate back solved x
  Dtilde[1] = Eigen::MatrixXd(nDof, NC);
  Eigen::MatrixXd N = Eigen::MatrixXd::Zero(nDof, NC), D[2];
  Eigen::VectorXd delta_i = Eigen::VectorXd::Zero(nDof);
  D[0] = Eigen::MatrixXd::Zero(nDof, NC); //concatenate back solved x
  D[1] = Eigen::MatrixXd::Zero(nDof, NC);
  std::vector < std::vector<Eigen::Vector3d> >x0s, v0s;
  for (size_t mi = 0; mi < meshes.size(); mi++){
    x0s.push_back(meshes[mi]->x);
    v0s.push_back(meshes[mi]->v);
  }

  int iters = 1;
  double residual = 0;
  double norm0 = 0;
  double xtol = 1e-7;
  int GS_ITER = 10 * std::min(10, NC);
  double fcscale = 4 / (h*h);
  double norm00 = -1;
  double norm_newton0 = -1;

  Eigen::VectorXd lambda = Eigen::VectorXd::Zero(NC);
  Eigen::VectorXd beta = Eigen::VectorXd::Zero(NC * 2);
  for (int mi = 0; mi < (int)meshes.size(); mi++){
    m = meshes[mi];
    updateStiffness();
    DampMats[mi] = DampMat;
    Ks[mi] = K;
    updateMb();
    Mbs[mi] = Mb;
  }
  
  //compute contact frames
  for (int ci = 0; ci < (int)collision.size(); ci++){
    Contact & c = collision[ci];
    ortho3(c.normal, c.d1, c.d2);
  }

  //vertex offset in global vector for each mesh.
  std::vector<int> voffset(meshes.size(), 0);
  for (int mi = 1; mi < (int)meshes.size(); mi++){
    voffset[mi] = (int)meshes[mi - 1]->x.size();
  }
  
  for (int mi = 0; mi < meshes.size(); mi++){
    m = meshes[mi];
    int offset = dim * voffset[mi];
    for (int ci = 0; ci < NC; ci++){
      Contact c = collision[ci];
      //solve for vertex response
      if (c.m[0] == mi){
        int vi = c.v1;
        for (int d = 0; d < dim; d++){
          N(offset + vi*dim + d, ci) = c.normal[d];
        }
        for (int d = 0; d < dim; d++){
          D[0](offset + vi*dim + d, ci) = c.d1[d];
        }
        for (int d = 0; d < dim; d++){
          D[1](offset + vi*dim + d, ci) = c.d2[d];
        }
      }
    }
  }

  bool newtonConverge = false;
  bool gsconverge = false;
  for (int newtonIter = 0; newtonIter < N_NEWTON; newtonIter++){
    for (int mi = 0; mi < meshes.size(); mi++){
      m = meshes[mi];
      int nDofi = dim * (int)m->x.size();
      x0 = x0s[mi];
      v0 = v0s[mi];
      delta = delta_i.block(dim * voffset[mi], 0, nDofi, 1);
      if (newtonIter > 0){
        updateStiffness();
        DampMats[mi] = DampMat;
        Ks[mi] = K;
      }
      else{
        DampMat = DampMats[mi];
        K = Ks[mi];
      }

      Mb = Mbs[mi];
      std::vector<double> x(nDofi);
      eigen2vector(delta, x);
      solver = (*solvers)[mi];
      //std::vector<double> rhs = getGrad();
      //Eigen::SparseMatrix<double> KM = getHessian();
      //setx(x.data());
      //solver->solve(KM, rhs, x, residual, iters);
      //vector2eigen(x, delta);
      newton->solver = solver;
      int maxIter0 = newton->maxIters;
      newton->maxIters = 1;
      int status = runNewton(x);
      if (status == -2){
        return status;
      }
      newton->maxIters = maxIter0;
      std::cout << "solved delta "<<newtonIter<<" " <<mi<<"\n";
      delta_i.block(dim * voffset[mi], 0, nDofi, 1) = delta;
      int offset = dim * voffset[mi];
      for (int ci = 0; ci < NC; ci++){
        Contact c = collision[ci];
        //solve for vertex response
        if (c.m[0] == mi){
          int vi = c.v1;
		  (*solvers)[mi]->solve(N.block(offset, ci, x.size(), 1).data(), x, residual, iters);
          Ntilde.block(offset, ci, x.size(), 1) = Eigen::Map<Eigen::VectorXd>(x.data(), x.size());

		  std::cout << "V:\n";
		  printVec(N.block(offset, ci, x.size(), 1).data(), x.size());
		  std::cout << "V tilde:\n";
		  printVec(x.data(), x.size());

		  (*solvers)[mi]->solve(D[0].block(offset, ci, x.size(), 1).data(), x, residual, iters);
          Dtilde[0].block(offset, ci, x.size(), 1) = Eigen::Map<Eigen::VectorXd>(x.data(), x.size());

		  std::cout << "D:\n";
		  printVec(D[0].block(offset, ci, x.size(), 1).data(), x.size());
		  std::cout << "D tilde:\n";
		  printVec(x.data(), x.size());

		  (*solvers)[mi]->solve(D[1].block(offset, ci, x.size(), 1).data(), x, residual, iters);
          Dtilde[1].block(offset, ci, x.size(), 1) = Eigen::Map<Eigen::VectorXd>(x.data(), x.size());
        }
      }
    }

    norm0 = -1;
    Eigen::VectorXd lambda0 = lambda, beta0 = beta;

    for (int gsIter = 0; gsIter < GS_ITER; gsIter++){
      for (int contacti = 0; contacti < lambda.rows(); contacti++){
        
        //pk = sum N_jl + D_jb +b j!=k
        Contact c = collision[contacti];
        int vi = c.v1;
		//vel in previous time step
		Eigen::Vector3d vt = v0[vi];
        Eigen::Vector3d pbar = delta.segment<3>(dim*vi);
        for (int j = 0; j < lambda.rows(); j++){
          if (j == contacti){
            continue;
          }
          for (int kk = 0; kk < dim; kk++){
            pbar[kk] += Ntilde(vi*dim + kk, j) * lambda(j);
            pbar[kk] += Dtilde[0](vi*dim + kk, j) * beta(j * 2);
            pbar[kk] += Dtilde[1](vi*dim + kk, j) * beta(j * 2 + 1);
          }
        }
        Eigen::Vector3d vbar = pbar;

        //no contact impulse.
        if (vbar.dot(c.normal) > 0){  //eq (14) > 0
          lambda[contacti] = 0;
          beta[2 * contacti] = 0;
          beta[2 * contacti + 1] = 0;
          continue;
        }
        Eigen::Vector3d ntilde = Eigen::Vector3d::Zero();
        Eigen::Vector3d dtilde1, dtilde2;
        for (int jj = 0; jj < dim; jj++){
          ntilde[jj] = Ntilde(vi*dim + jj, contacti);
          dtilde1[jj] = Dtilde[0](vi*dim + jj, contacti);
          dtilde2[jj] = Dtilde[1](vi*dim + jj, contacti);
        }
        Eigen::Matrix3d A;

        Eigen::Vector3d n, d1, d2, b;
        n = c.normal;
        d1 =c.d1;
        d2 =c.d2;
        A << n.dot(ntilde), n.dot(dtilde1), n.dot(dtilde2),
          d1.dot(ntilde), d1.dot(dtilde1), d1.dot(dtilde2),
          d2.dot(ntilde), d2.dot(dtilde1), d2.dot(dtilde2);
        b << n.dot(vbar), d1.dot(vbar), d2.dot(vbar);

        //used for other experiments
        double mu = 0.3; // TODO: make mu a config file term passed in;
        //walker;
        //mu = 1;
        // TODO: change to mu[contacti] to handle e.g., treat plastic/plastic and plastic/metal

        // 1. test break condition:
        // done above

        // 2. test stick condition:
        Eigen::Vector3d x = -A.inverse()*b;
        lambda[contacti] = x[0];
        beta[2 * contacti] = x[1];
        beta[2 * contacti + 1] = x[2];
        double coulombCondition = x[0] * x[0] * mu * mu >= x[1] * x[1] + x[2] * x[2];
        if ((frictionModel == FINITE) && !coulombCondition){
          // 3. Else must be slide condition: solve w fixed point.
          double lambda_i = lambda[contacti];
          Eigen::Vector2d beta_i = Eigen::Vector2d::Zero();
          Eigen::Vector2d betaOld_i = beta_i;
          int max_iter = 100;
          int iter = 0;
          for (iter = 0; iter < max_iter; iter++){
            betaOld_i = beta_i;
            // lambda_i = -n.dot(vbar + Tcomp * beta_i) / n.dot(ncomp);
            Eigen::Vector3d rf = dtilde1 * beta_i[0] + dtilde2 * beta_i[1];
            lambda_i = -(b(0) + n.dot(rf)) / A(0, 0);

			auto Afriction = A.block(1, 1, 2, 2);
			Eigen::Vector2d bfriction(b[1],b[2]);
			Eigen::Vector3d r = ntilde* lambda_i - 0.5*h*vt;
			bfriction[0] += d1.dot(r);
			bfriction[1] += d2.dot(r);
			beta_i = -0.5 * Afriction.inverse()*bfriction;
            // beta_i = -T.transpose() * (vbar + Tcomp *beta_i + ncomp * lambda_i);
          
            beta_i.normalize();
            beta_i *= mu * lambda_i;
            //std::cout << betaOld_i << " " << beta_i << " beta\n";
            if ((betaOld_i - beta_i).norm() < 1e-12){
              break;
            }
          }

         lambda[contacti] = lambda_i;
         beta[2 * contacti] = beta_i(0);
         beta[2 * contacti + 1] = beta_i(1);
        }
      }
      double norm1 = 0;
	  std::cout << "beta:\n";
	  for (int ci = 0; ci < NC; ci++){
        norm1 = std::max(norm1, fcscale*std::abs(lambda[ci] - lambda0[ci]));
        norm1 = std::max(norm1, fcscale*std::abs(beta[2 * ci] - beta0[2 * ci]));
        norm1 = std::max(norm1, fcscale*std::abs(beta[2 * ci + 1] - beta0[2 * ci + 1]));
		std::cout << beta[2 * ci] << " ";
		std::cout << beta[2 * ci+1] << " \n";
      }

      //std::cout << "newton gs iter " << newtonIter << " " << gsIter << "====\n";
      //std::cout << norm0 << " " << norm1 << "\n";
      beta0 = beta;
      lambda0 = lambda;
      norm0 = norm1;
      if (gsIter == 0){
        norm00 = norm1;
      }
      if (gsIter>5 && norm1 < 1e-8){
        gsconverge = true;
        break;
      }
    }

    for (int ci = 0; ci < NC; ci++){
      for (int ii = 0; ii < nDof; ii++){
        delta_i[ii] += lambda[ci] * Ntilde(ii, ci);
        delta_i[ii] += beta[2 * ci] * Dtilde[0](ii, ci);
        delta_i[ii] += beta[2 * ci + 1] * Dtilde[1](ii, ci);
      }
    }
    //maximum change in delta.
    //used for convergence check.
    double maxDiff = 0;
    //check applied impulse
    std::vector<Eigen::Vector3d> totalImpulse(meshes.size(), Eigen::Vector3d::Zero());
    for (int mi = 0; mi < (int)meshes.size(); mi++){
      m = meshes[mi];
      int nDofi = dim * (int)m->x.size();
      for (unsigned int ii = 0; ii < m->X.size(); ii++){
        m->fc[ii] = Eigen::Vector3d::Zero();
      }
      delta = delta_i.block(dim * voffset[mi], 0, nDofi, 1);
      //update contact force for rendering
      for (int ci = 0; ci < NC; ci++){
        for (int ii = 0; ii < nDofi; ii++){
          int vi = ii / 3;
          //m->fc[vi][ii % 3] += fcscale * lambda[ci] * Ntilde(dim * voffset[mi] + ii, ci);
          //m->fc[vi][ii % 3] += fcscale * beta[2 * ci] * Dtilde[0](dim * voffset[mi] + ii, ci);
          //m->fc[vi][ii % 3] += fcscale * beta[2 * ci + 1] * Dtilde[1](dim * voffset[mi] + ii, ci);
          m->fc[vi][ii % 3] += fcscale * lambda[ci] * N(dim * voffset[mi] + ii, ci);
          m->fc[vi][ii % 3] += fcscale * beta[2 * ci] * D[0](dim * voffset[mi] + ii, ci);
          m->fc[vi][ii % 3] += fcscale * beta[2 * ci + 1] * D[1](dim * voffset[mi] + ii, ci);
        }
      }
      x0 = x0s[mi];
      v0 = v0s[mi];
      //update position and velocity using delta
      for (unsigned int ii = 0; ii < m->x.size(); ii++){
        for (int jj = 0; jj < dim; jj++){
          m->x[ii][jj] = x0[ii][jj] + delta[ii*dim + jj];
        }
        for (int jj = 0; jj < dim; jj++){
          m->v[ii][jj] = (2.0 / h) * delta[ii*dim + jj] - v0[ii][jj];
        }
      }
      Mb = Mbs[mi];
      DampMat = DampMats[mi];
      std::vector<double> grad = getGrad();
      for (int ci = 0; ci < NC; ci++){
        for (int vi = 0; vi < (int)m->x.size(); vi++){
          for (int di = 0; di < dim; di++){
            grad[dim * vi + di] += lambda[ci] * N(dim * (voffset[mi] + vi) + di, ci);
            grad[dim * vi + di] += beta[2 * ci] * D[0](dim * (voffset[mi] + vi) + di, ci);
            grad[dim * vi + di] += beta[2 * ci + 1] * D[1](dim * (voffset[mi] + vi) + di, ci);
            //totalImpulse[mi][di] += lambda[ci] * N(dim * (voffset[mi] + vi) + di, ci);
            //totalImpulse[mi][di] += beta[2 * ci] * D[0](dim * (voffset[mi] + vi) + di, ci);
            //totalImpulse[mi][di] += beta[2 * ci + 1] * D[1](dim * (voffset[mi] + vi) + di, ci);
          }
        }
      }
      for (int ii = 0; ii < (int)grad.size(); ii++){
        maxDiff = std::max(maxDiff, std::abs(grad[ii]));
      }
      //std::cout << "totalImpulse " << mi << " " << totalImpulse[mi][0] << " " << totalImpulse[mi][1] << " " << totalImpulse[mi][2] << "\n";
    }
    //std::cout << "newton " << newtonIter << " " << norm_newton0 << " " << maxDiff << " " << maxDiff / norm_newton0 << "\n";
    if (newtonIter > 0){
      if (maxDiff / norm_newton0 < 1e-9 || maxDiff < 1e-19){
        break;
      }
    }
    else{
      norm_newton0 = maxDiff;
    }
  }
  //if (!gsconverge){
  //  int input;
  //  std::cin >> input;
  //}
  //remove normal component of contact vertex
  //double E0 = m->getAllEnergy(true);
  //only works for 1 mesh currently.
  Eigen::VectorXd y = Eigen::VectorXd::Zero(collision.size());
  Eigen::VectorXd b(collision.size());
  Eigen::MatrixXd A(collision.size(), collision.size());
  for (int i = 0; i < b.rows(); i++){
    int vi = collision[i].v1;
    Eigen::Vector3d N = collision[i].normal;
    b[i] = m->v[vi].transpose() * N;
    for (int j = 0; j < b.rows(); j++){
      int vj = collision[j].v1;
      Eigen::Vector3d ntilde = Ntilde.block(dim * vj, i, 3, 1);
      A(i, j) = N.transpose() * ntilde;
    }
  }
  //std::cout << A << "\n";
  //std::cout << b << "\n";
  Eigen::VectorXd y0 = y;
  bool velProjection = false;
  //bool velProjection = true;
  if (velProjection){
    int projectionIter = std::min(20, (int)y.rows());
    double projectionThresh = 1e-15;
    for (int iter = 0; iter < projectionIter; iter++){
      y0 = y;
      for (int i = 0; i < collision.size(); i++){
        double bi = b[i];
        for (int j = 0; j < y.rows(); j++){
          if (j == i){
            continue;
          }
          bi -= A(i, j) * y(j);
        }
        y(i) = bi / A(i, i);
        //if (y(i) < 0){
        //  y(i) = 0;
        //}
      }
      double norm = (y0 - y).norm() / y.rows();
      //std::cout << iter << " projection gs norm " << norm << " .\n";
      if (norm < projectionThresh){
        break;
      }
    }
    //y = A.ldlt().solve(b);
    for (int i = 0; i < collision.size(); i++){
      //std::cout << y(i) << "\n";
      for (size_t j = 0; j < m->v.size(); j++){
        for (int k = 0; k < dim; k++){
          m->v[j][k] -= y(i) * Ntilde(dim * j + k, i);
        }
      }
    }
  }
  //for (int i = 0; i < collision.size(); i++){
  //  std::cout << m->v[collision[i].v1][1] << "\n";
  //}
//old projection code.
  bool oldProjection = false;
  //bool oldProjection = true;
  if (oldProjection){
    m->vu = m->v;
    for (int i = 0; i<collision.size(); i++){
      int vi = collision[i].v1;
      if (lambda[i]>1e-20){
        m->v[vi] = m->v[vi] - m->v[vi].dot(collision[i].normal) * collision[i].normal;
      }
    }
  }
  //double E1 = m->getAllEnergy(true);
  //std::cout << " Energies " << E0 << " " << E1 << "\n";
  return 0;
}

///\param ci contact index.
///\param lambda output array of lambda.
void solveFriction(Stepper::FrictionModel frictionModel,
  int ci, const Contact & c,
  Eigen::Vector3d vbar, double mu,
  std::vector<Eigen::SparseVector<double> > & Nsp,
  std::vector<Eigen::SparseVector<double> > * Dsp,
  const Eigen::MatrixXd & Ntilde, const Eigen::MatrixXd * Dtilde,
  Stepper * stepper,
  double * lambda, double *beta)
{
  float h = stepper->h;
  Eigen::Matrix3d A;
  Eigen::Vector3d b;
  A << Nsp[ci].dot(Ntilde.col(ci)), Nsp[ci].dot(Dtilde[0].col(ci)), Nsp[ci].dot(Dtilde[1].col(ci)),
    Dsp[0][ci].dot(Ntilde.col(ci)), Dsp[0][ci].col(ci).dot(Dtilde[0].col(ci)), Dsp[0][ci].dot(Dtilde[1].col(ci)),
    Dsp[1][ci].dot(Ntilde.col(ci)), Dsp[1][ci].col(ci).dot(Dtilde[0].col(ci)), Dsp[1][ci].col(ci).dot(Dtilde[1].col(ci));
  b << c.normal.dot(vbar), c.d1.dot(vbar), c.d2.dot(vbar);

  //add frictionless solve
  //lambda = N dot ntilde/vbar dot ntilde
  if (frictionModel == Stepper::INFINITE){
    Eigen::Vector3d x = -A.inverse()*b;
    lambda[ci] = x[0];
    beta[2 * ci] = x[1];
    beta[2 * ci+ 1] = x[2];
  }
  else if (frictionModel == Stepper::NONE){
    lambda[ci] = -b(0, 0) / A(0, 0);
    beta[2 * ci] = 0;
    beta[2 * ci + 1] = 0;
  }
  else if (frictionModel == Stepper::TRESCA){
    double mu_tresca = 1e-6; // tresca friction proportion
    Eigen::Vector3d x = -A.inverse()*b;
    //Eigen::Vector3d x;
    x[0] = -b(0, 0) / A(0, 0);
    lambda[ci] = x[0];
    beta[2 * ci] = 0;
    beta[2 * ci + 1] = 0;
    Eigen::Vector2d f = -A.block(1, 1, 2, 2).inverse()*b.segment<2>(1);
    if (f.norm() > 1e-9){
      f.normalize();
      f *= mu_tresca * h * h;
      beta[2 * ci] = f[0];
      beta[2 * ci + 1] = f[1];
    }
  }
  else if (frictionModel == Stepper::FINITE){
    // 1. test break condition:
    // done above

    // 2. test stick condition:
    Eigen::Vector3d x = -A.inverse()*b;
    lambda[ci] = x[0];
    beta[2 * ci] = x[1];
    beta[2 * ci + 1] = x[2];
    double coulombCondition = x[0] * x[0] * mu * mu >= x[1] * x[1] + x[2] * x[2];

    if (!coulombCondition){
      // 3. Else must be slide condition: solve w fixed point.

      double lambda_i = 0;
      Eigen::Vector2d beta_i = Eigen::Vector2d::Zero();
      Eigen::Vector2d betaOld_i = beta_i;
      int max_iter = 30;
      int iter = 0;
      for (iter = 0; iter < max_iter; iter++){
        betaOld_i = beta_i;
        // lambda_i = -n.dot(vbar + Tcomp * beta_i) / n.dot(ncomp);
        Eigen::VectorXd rf = Dtilde[0].col(ci) * beta_i(0) + Dtilde[1].col(ci) * beta_i(1);
        lambda_i = -(b(0) + Nsp[ci].dot(rf)) / A(0, 0);

        // beta_i = -T.transpose() * (vbar + Tcomp *beta_i + ncomp * lambda_i);
        Eigen::VectorXd r = rf + Ntilde.col(ci) * lambda_i;
        beta_i(0) = -b(1) - Dsp[0][ci].dot(r);
        beta_i(1) = -b(2) - Dsp[1][ci].dot(r);

        beta_i = mu * lambda_i * beta_i.normalized();
        if ((betaOld_i - beta_i).norm() < 1e-12){
          break;
        }
      }

      lambda[ci] = lambda_i;
      beta[2 * ci] = beta_i(0);
      beta[2 * ci + 1] = beta_i(1);
    }
  }
}

double dot(Eigen::SparseVector<double> & s, Eigen::VectorXd v){
  double sum = 0;
  Eigen::SparseVector<double>::InnerIterator it(s);
  for (; it; ++it){
    sum += it.value() * v[it.index()];
  }
  return sum;
}

void
StepperNewmark::resolveCollision3D(std::vector<Contact> & collision)
{
  if (useOldContact){
    resolveCollision3D_old(collision);
    return;
  }
  int dim = m->dim;
  int N_NEWTON = 4;
  std::vector<Eigen::SparseMatrix<double> >DampMats(meshes.size()), Ks(meshes.size());
  std::vector<std::vector<Eigen::Vector3d> > Mbs(meshes.size());
  //std::cout << collision.size() << " contacts.\n";
  //number of constraints
  int NC = (int)collision.size();
  //total dof
  int nDof = 0;
  for (int mi = 0; mi < meshes.size(); mi++){
    m = meshes[mi];
    nDof += dim * (int)m->X.size();
  }
  //back solve to compute Ntilde and Dtilde
  Eigen::MatrixXd Ntilde(nDof, NC), Dtilde[2];
  Dtilde[0] = Eigen::MatrixXd(nDof, NC);
  Dtilde[1] = Eigen::MatrixXd(nDof, NC);
  Eigen::MatrixXd N = Eigen::MatrixXd::Zero(nDof, NC), D[2];
  std::vector<Eigen::SparseVector<double> > Nsp (NC), Dsp[2];
  Eigen::VectorXd delta_i = Eigen::VectorXd::Zero(nDof);
  D[0] = Eigen::MatrixXd::Zero(nDof, NC);
  D[1] = Eigen::MatrixXd::Zero(nDof, NC);
  Dsp[0] = std::vector<Eigen::SparseVector<double> > (NC);
  Dsp[1] = std::vector<Eigen::SparseVector<double> > (NC);

  std::vector < std::vector<Eigen::Vector3d> >x0s, v0s;
  for (size_t mi = 0; mi < meshes.size(); mi++){
    x0s.push_back(meshes[mi]->x);
    v0s.push_back(meshes[mi]->v);
  }

  int iters = 1;
  double residual = 0;
  double norm0 = 0;
  double xtol = 1e-7;
  int GS_ITER = 20 * std::min(10, NC);
  double fcscale = 4 / (h*h);
  double norm00 = -1;
  double norm_newton0 = -1;
  Eigen::VectorXd lambda = Eigen::VectorXd::Zero(NC);
  Eigen::VectorXd beta = Eigen::VectorXd::Zero(NC * 2);
  for (int mi = 0; mi < (int)meshes.size(); mi++){
    m = meshes[mi];
    updateStiffness();
    DampMats[mi] = DampMat;
    Ks[mi] = K;
    updateMb();
    Mbs[mi] = Mb;
  }
  //compute contact frames
  for (int ci = 0; ci < (int)collision.size(); ci++){
    Contact & c = collision[ci];
    ortho3(c.normal, c.d1, c.d2);
  }
  //vertex offset in global vector for each mesh.
  std::vector<int> voffset(meshes.size(),0);
  for (int mi = 1; mi < (int)meshes.size(); mi++){
    voffset[mi] = (int)meshes[mi - 1]->x.size();
  }

  //compute constant normal and tangent contact vectors.
  for (int ci = 0; ci < NC; ci++){
    Nsp[ci] = Eigen::SparseVector<double>(nDof);
    Dsp[0][ci] = Eigen::SparseVector<double>(nDof);
    Dsp[1][ci] = Eigen::SparseVector<double>(nDof);
  }
  for (int mi = 0; mi < meshes.size(); mi++){
    int offset = dim * voffset[mi];
    for (int ci = 0; ci < NC; ci++){
      Contact c = collision[ci];
      //vertex response
      if (c.m[0] == mi){
        int vi = c.v1;
        for (int d = 0; d < dim; d++){
          int dof = offset + vi*dim + d;
          N(dof, ci) = c.normal[d];
          Nsp[ci].insert(dof) = c.normal[d];
        }
        for (int d = 0; d < dim; d++){
          int dof = offset + vi*dim + d;
          D[0](dof, ci) = c.d1[d];
          Dsp[0][ci].insert(dof) = c.d1[d];
        }
        for (int d = 0; d < dim; d++){
          int dof = offset + vi*dim + d;
          D[1](dof, ci) = c.d2[d];
          Dsp[1][ci].insert(dof) = c.d2[d];
        }
      }
      //quad response
      else if (c.m[1] == mi){
        float W[4];
        bilinearWeights(c.alpha, W);
        //loop over 4 vertices of a quad
        for (int j = 0; j < 4; j++){
          int vi = c.v2[j];
          Eigen::Vector3d normal = W[j] * c.normal;
          for (int d = 0; d < dim; d++){
            int dof = offset + vi*dim + d;
            N(dof, ci) = -normal[d];
            Nsp[ci].insert(dof) = -normal[d];
          }
        }
        for (int j = 0; j < 4; j++){
          int vi = c.v2[j];
          for (int d = 0; d < dim; d++){
            int dof = offset + vi*dim + d;
            D[0](dof, ci) = -c.d1[d];
            Dsp[0][ci].insert(dof) = -c.d1[d];
          }
        }
        for (int j = 0; j < 4; j++){
          int vi = c.v2[j];
          for (int d = 0; d < dim; d++){
            int dof = offset + vi*dim + d;
            D[1](dof, ci) = -c.d2[d];
            Dsp[1][ci].insert(dof) = -c.d2[d];
          }
        }
      }
    }
  }

  bool newtonConverge = false;
  bool gsconverge = false;
  for (int newtonIter = 0; newtonIter < N_NEWTON; newtonIter++){
    for (int mi = 0; mi < meshes.size(); mi++){
      m = meshes[mi];
      int nDofi = dim * (int)m->x.size();
      x0 = x0s[mi];
      v0 = v0s[mi];
      delta = delta_i.block(dim * voffset[mi], 0, nDofi, 1);
      if (newtonIter > 0){
        updateStiffness();
        DampMats[mi] = DampMat;
        Ks[mi] = K;
      }
      else{
        DampMat = DampMats[mi];
        K = Ks[mi];
      }
      
      Mb = Mbs[mi];
      std::vector<double> x(nDofi);
      eigen2vector(delta, x);
      solver = (*solvers)[mi];
      //std::vector<double> rhs = getGrad();
      //Eigen::SparseMatrix<double> KM = getHessian();
      //setx(x.data());
      //solver->solve(KM, rhs, x, residual, iters);
      //vector2eigen(x, delta);
      newton->solver = solver;
      int maxIter0 = newton->maxIters;
      newton->maxIters = 1;
      int status = runNewton(x);
      newton->maxIters = maxIter0;

      delta_i.block(dim * voffset[mi], 0, nDofi, 1) = delta;
      int offset = dim * voffset[mi];
      //solve for compliant contact response
      for (int ci = 0; ci < NC; ci++){
          (*solvers)[mi]->solve(N.block(offset, ci, x.size(), 1).data(), x, residual, iters);
          Ntilde.block(offset, ci, x.size(), 1) = Eigen::Map<Eigen::VectorXd>(&(x[0]), x.size());
          (*solvers)[mi]->solve(D[0].block(offset, ci, x.size(), 1).data(), x, residual, iters);
          Dtilde[0].block(offset, ci, x.size(), 1) = Eigen::Map<Eigen::VectorXd>(&(x[0]), x.size());
          (*solvers)[mi]->solve(D[1].block(offset, ci, x.size(), 1).data(), x, residual, iters);
          Dtilde[1].block(offset, ci, x.size(), 1) = Eigen::Map<Eigen::VectorXd>(&(x[0]), x.size());
      }
    }

    norm0 = -1;
    Eigen::VectorXd lambda0 = lambda, beta0 = beta;

    for (int gsIter = 0; gsIter < GS_ITER; gsIter++){
      for (int contacti = 0; contacti < lambda.rows(); contacti++){
        //pk = sum N_jl + D_jb +b j!=k
        Contact c = collision[contacti];
        Eigen::VectorXd pbar = delta_i;
        for (int j = 0; j < lambda.rows(); j++){
          if (j == contacti){
            continue;
          }
          pbar += Ntilde.col(j) * lambda(j);
          pbar += Dtilde[0].col(j) * beta(j * 2);
          pbar += Dtilde[1].col(j) * beta(j * 2 + 1);
        }
        double vn = dot(Nsp[contacti],pbar);
        double vt1 = dot(Dsp[0][contacti], pbar);
        double vt2 = dot(Dsp[1][contacti], pbar);

        Eigen::Vector3d vbar = vn * c.normal + vt1 * c.d1 + vt2*c.d2;

        //no contact impulse.
        if (vn>0){  //eq (14) > 0
          lambda[contacti] = 0;
          beta[2 * contacti] = 0;
          beta[2 * contacti + 1] = 0;
          continue;
        }
        double mu = 0.7; // TODO: make mu a config file term passed in;
        // TODO: change to mu[contacti] to handle e.g., treat plastic/plastic and plastic/metal

        solveFriction(frictionModel, contacti, c, vbar, mu,
          Nsp,Dsp, Ntilde, Dtilde, this, lambda.data(), beta.data());
			}
      double norm1 = 0;
      for (int ci = 0; ci < NC; ci++){
        norm1 = std::max(norm1, fcscale*std::abs(lambda[ci] - lambda0[ci]));
        norm1 = std::max(norm1, fcscale*std::abs(beta[2 * ci] - beta0[2 * ci]));
        norm1 = std::max(norm1, fcscale*std::abs(beta[2 * ci + 1] - beta0[2 * ci + 1]));
      }
      //std::cout<<"newton gs iter "<<newtonIter<<" "<<gsIter<<"====\n";
      //std::cout<<norm0<<" "<<norm1<<"\n";
      beta0 = beta;
      lambda0 = lambda;
      norm0 = norm1;
      if (gsIter == 0){
        norm00 = norm1;
      }
      if (gsIter>5 && norm1 < 1e-8){//gsIter>5 &&
        gsconverge = true;
        break;
      }
    }

    for (int ci = 0; ci < NC; ci++){
      for (int ii = 0; ii < nDof; ii++){
        delta_i[ii] += lambda[ci] * Ntilde(ii, ci);
        delta_i[ii] += beta[2 * ci] * Dtilde[0](ii, ci);
        delta_i[ii] += beta[2 * ci + 1] * Dtilde[1](ii, ci);
      }
    }
    //maximum change in delta.
    //used for convergence check.
    double maxDiff = 0;
    //check applied impulse
    std::vector<Eigen::Vector3d> totalImpulse(meshes.size(), Eigen::Vector3d::Zero());
    for (int mi = 0; mi < (int)meshes.size(); mi++){
      m = meshes[mi];
      int nDofi = dim * (int)m->x.size();
      for (unsigned int ii = 0; ii < m->X.size(); ii++){
        m->fc[ii] = Eigen::Vector3d::Zero();
      }
      delta = delta_i.block(dim * voffset[mi], 0, nDofi, 1);
      //update contact force for rendering
      for (int ci = 0; ci < NC; ci++){
        for (int ii = 0; ii < nDofi; ii++){
          int vi = ii / 3;
          m->fc[vi][ii % 3] += fcscale * lambda[ci] * Ntilde(dim * voffset[mi] + ii, ci);
          m->fc[vi][ii % 3] += fcscale * beta[2 * ci] * Dtilde[0](dim * voffset[mi] + ii, ci);
          m->fc[vi][ii % 3] += fcscale * beta[2 * ci + 1] * Dtilde[1](dim * voffset[mi] + ii, ci);
          //m->fc[vi][ii % 3] += fcscale * lambda[ci] * N(dim * voffset[mi] + ii, ci);
          //m->fc[vi][ii % 3] += fcscale * beta[2 * ci] * D[0](dim * voffset[mi] + ii, ci);
          //m->fc[vi][ii % 3] += fcscale * beta[2 * ci + 1] * D[1](dim * voffset[mi] + ii, ci);
        }
      }
      x0 = x0s[mi];
      v0 = v0s[mi];
      //update position and velocity using delta
      for (unsigned int ii = 0; ii < m->x.size(); ii++){
        for (int jj = 0; jj < dim; jj++){
          m->x[ii][jj] = x0[ii][jj] + delta[ii*dim + jj];
        }
        for (int jj = 0; jj < dim; jj++){
          m->v[ii][jj] = (2.0 / h) * delta[ii*dim + jj] - v0[ii][jj];
        }
      }
      Mb = Mbs[mi];
      DampMat = DampMats[mi];
      std::vector<double> grad = getGrad();
      for (int ci = 0; ci < NC; ci++){
        for (int vi = 0; vi < (int)m->x.size(); vi++){
          for (int di = 0; di < dim; di++){
            grad[dim * vi + di] += lambda[ci] * N(dim * (voffset[mi] + vi) + di, ci);
            grad[dim * vi + di] += beta[2*ci] * D[0](dim * (voffset[mi] + vi) + di, ci);
            grad[dim * vi + di] += beta[2 * ci + 1] * D[1](dim * (voffset[mi] + vi) + di, ci);
            //totalImpulse[mi][di] += lambda[ci] * N(dim * (voffset[mi] + vi) + di, ci);
            //totalImpulse[mi][di] += beta[2 * ci] * D[0](dim * (voffset[mi] + vi) + di, ci);
            //totalImpulse[mi][di] += beta[2 * ci + 1] * D[1](dim * (voffset[mi] + vi) + di, ci);
          }
        }
      }
      for (int ii = 0; ii <(int)grad.size(); ii++){
        maxDiff = std::max(maxDiff, std::abs(grad[ii]));
      }
      //std::cout << "totalImpulse " << mi << " " << totalImpulse[mi][0] << " " << totalImpulse[mi][1] << " " << totalImpulse[mi][2] << "\n";
    }
    //std::cout << "newton " << newtonIter << " " << norm_newton0 << " " << maxDiff << " " << maxDiff / norm_newton0 << "\n";
    if (newtonIter > 0){
      if (maxDiff / norm_newton0 < 1e-9 || maxDiff < 1e-19){
        break;
      }
    }
    else{
      norm_newton0 = maxDiff;
    }
  }
  //if (!gsconverge){
  //  int input;
  //  std::cin >> input;
  //}

}
