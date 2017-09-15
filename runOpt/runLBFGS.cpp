#include "runLBFGS.hpp"
#include <stdio.h>

///numerical gradient
Eigen::VectorXd numGrad(float h, RealFun & f,
  Eigen::VectorXd x)
{
  std::cout << "Num grad ======================\n";
  Eigen::VectorXd grad(x.size());
  for (unsigned int ii = 0; ii<x.size(); ii++){
    x[ii] += h;
    f.setParam(x);
    double fp = f.f();

    x[ii] -= 2 * h;
    f.setParam(x);
    double fm = f.f();

    x[ii] += h;
    std::cout << fm << " " << fp << "\n";
    grad[ii] = (fp - fm) / (2 * h);
  }
  std::cout << "Num grad end ======================\n";
  return grad;
}


int runLBFGS(const ConfigFile & config,
  JumperObjective &objective)
{
  //std::cout<<"Optimize design parameters.\n";
  //int N = 2;
  //Eigen::VectorXd param(N);
  //for(int ii = 0; ii<param.rows(); ii++){
  //  param[ii] = 0;
  //}

  //objective.init(param);
  //objective.setParam(param);

  //lbfgsfloatval_t fx;
  //lbfgsfloatval_t *m_x = lbfgs_malloc(N);
  //lbfgs_parameter_t lbfgs_param;
  //lbfgs_parameter_init(&lbfgs_param);
  //lbfgs_param.max_step = 5e-6;
  //if (m_x == NULL) {
  //    printf("ERROR: Failed to allocate a memory block for variables.\n");
  //    return 1;
  //}

  ///* Initialize the variables. */
  //for (int i = 0;i < N;i ++) {
  //    m_x[i] = 0;
  //}

  int ret = 0;
  //ret = lbfgs(N, m_x, &fx, _evaluate, _progress, &objective, &lbfgs_param);
  //printf("L-BFGS optimization terminated with status code = %d\n", ret);
  //printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, m_x[0], m_x[1]);

  return ret;
}



lbfgsfloatval_t _evaluate(
  void *instance,
  const lbfgsfloatval_t *x,
  lbfgsfloatval_t *g,
  const int n,
  const lbfgsfloatval_t step
  )
{
  //    lbfgsfloatval_t fx = 0.0;

  //    lbfgsfloatval_t t1 = 2*(1 - x[0]);
  //    lbfgsfloatval_t t2 = 2 - x[1];
  //    lbfgsfloatval_t t3 = (x[1] - x[0]);
  //    fx += t1 * t1 + t2 * t2 + t3 * t3;

  //    g[0] = -2*t1 -2*t3;
  //    g[1] = -2*t2 + 2*t3;

  RealFun * obj = (RealFun*)(instance);
  double h = 0.0001*0.0254;
  Eigen::VectorXd param(n);
  for (int ii = 0; ii<n; ii++){
    param[ii] = x[ii];
  }
  Eigen::VectorXd grad = numGrad(h, *obj, param);
  if (n >= 2){
    std::cout << "iteration: " << step << " grad " << grad[0] << " " << grad[1] << "\n";
    std::cout << "Param " << param[0] << " " << param[1] << "\n";
  }

  obj->setParam(param);
  lbfgsfloatval_t val = obj->f();
  for (int ii = 0; ii<grad.size(); ii++){
    g[ii] = grad[ii];
  }
  return val;
}


int _progress(
  void *instance,
  const lbfgsfloatval_t *x,
  const lbfgsfloatval_t *g,
  const lbfgsfloatval_t fx,
  const lbfgsfloatval_t xnorm,
  const lbfgsfloatval_t gnorm,
  const lbfgsfloatval_t step,
  int n,
  int k,
  int ls
  )
{
  if (n >= 2){
    printf("Iteration %d:\n", k);
    printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    printf("\n");
  }
  return 0;
}

