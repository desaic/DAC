// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id: MyNLP.cpp 2005 2011-06-06 12:55:16Z stefan $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-11-05

#include "IpInterface.hpp"
#include "RealFun.hpp"
#include <iostream>
#include <sstream>
#include <cassert>

/* Ipopt Function Declarations */
Bool eval_fj(Index n, Number* x, Bool new_x,
  Number* obj_value, UserDataPtr user_data);

Bool eval_grad_fj(Index n, Number* x, Bool new_x,
  Number* grad_f, UserDataPtr user_data);

Bool eval_gj(Index n, Number* x, Bool new_x,
  Index m, Number* g, UserDataPtr user_data);

Bool eval_jac_gj(Index n, Number *x, Bool new_x,
  Index m, Index nele_jac,
  Index *iRow, Index *jCol, Number *values,
  UserDataPtr user_data);

Bool eval_hj(Index n, Number *x, Bool new_x, Number obj_factor,
  Index m, Number *lambda, Bool new_lambda,
  Index nele_hess, Index *iRow, Index *jCol,
  Number *values, UserDataPtr user_data);

Bool intermediate_cbj(Index alg_mod, Index iter_count, Number obj_value,
  Number inf_pr, Number inf_du, Number mu, Number d_norm,
  Number regularization_size, Number alpha_du,
  Number alpha_pr, Index ls_trials, UserDataPtr user_data);

/* Constructor. */
IpInterface::IpInterface() :objective(0),
nlp(0), logger(0)
{}

IpInterface::~IpInterface()
{}

int IpInterface::init()
{
  Index n = (Index)objective->param.size();
  //no Ipopt constraints
  Index m = 0;
  Index index_style = 0;
  if (x_L.size() != n){
    std::cout << "IpInterface un-init x_L\n";
    x_L.resize(n, -1);
  }
  if (x_U.size() != n){
    std::cout << "IpInterface un-init x_L\n";
    x_U.resize(n, 1);
  }
  nlp = CreateIpoptProblem(n, x_L.data(), x_U.data(), m, NULL, NULL, 0, 0,
    index_style, &eval_fj, &eval_gj, &eval_grad_fj,
    &eval_jac_gj, &eval_hj);
  AddIpoptNumOption(nlp, "tol", 1e-4);
  AddIpoptStrOption(nlp, "mu_strategy", "adaptive");
  AddIpoptStrOption(nlp, "hessian_approximation", "limited-memory");
  AddIpoptIntOption(nlp, "accept_after_max_steps", 10);
  //AddIpoptIntOption(nlp, "max_soc", 0);
  //AddIpoptIntOption(nlp, "watchdog_shortened_iter_trigger", 0);
  //AddIpoptStrOption(nlp, "derivative_test", "first-order");
  //AddIpoptNumOption(nlp, "derivative_test_perturbation", 2e-5);
  SetIntermediateCallback(nlp, intermediate_cbj);
  has_f_cache = false;
  has_grad_cache = false;
  p_cache=Eigen::VectorXd::Zero(n);
  return 0;
}

int IpInterface::run()
{
  Number objVal;
  Index i;
  Index n = (Index)objective->param.size();
  std::vector<Number>x(objective->param.size(), 0);
  //copy initial point from objective.
  for (size_t i = 0; i < x.size(); i++){
    x[i] = objective->param[i];
  }
  std::vector<Number>mult_xl(objective->param.size(), 0);
  std::vector<Number>mult_xu(objective->param.size(), 0);
  ApplicationReturnStatus status = IpoptSolve(nlp, x.data(), NULL, &objVal, NULL, 
    mult_xl.data(), mult_xu.data(), this);
  if (status == Solve_Succeeded) {
    printf("\n\nSolution of the primal variables, x\n");
    for (i = 0; i<n; i++) {
      printf("x[%d] = %e\n", i, x[i]);
    }

    printf("\n\nSolution of the bound multipliers, z_L and z_U\n");
    for (i = 0; i<n; i++) {
      printf("z_L[%d] = %e\n", i, mult_xl[i]);
    }
    for (i = 0; i<n; i++) {
      printf("z_U[%d] = %e\n", i, mult_xu[i]);
    }

    printf("\n\nObjective value\n");
    printf("f(x*) = %e\n", objVal);
  }
  else {
    printf("\n\nERROR OCCURRED DURING IPOPT OPTIMIZATION.\n");
  }
  return 0;
}


/* Function Implementations */
Bool eval_fj(Index n, Number* x, Bool new_x,
  Number* obj_value, UserDataPtr user_data)
{
  IpInterface * j = (IpInterface*)user_data;
  std::ostringstream ss;
  ss << "eval ";
  j->p_cache.resize(n);
  Eigen::VectorXd p(n);
  for (int ii = 0; ii<n; ii++){
    p[ii] = x[ii];
    ss << x[ii] << " ";
  }
  Eigen::VectorXd diff = j->p_cache - p;
  j->p_cache = p;
  if (diff.squaredNorm()<1e-12){
    new_x = false;
  }
  if (!new_x && j->has_f_cache){
    *obj_value = (Number)j->f_cache;
    ss << "\nval " << *obj_value << "\n";
    j->logMsg(ss.str());
    return true;
  }
  j->logMsg(ss.str());
  ss.clear();

  j->objective->setParam(p);
  *obj_value = j->objective->f();
    j->f_cache = *obj_value;
  j->has_f_cache = true;
  if (new_x){
    j->has_grad_cache = false;
  }

  ss << "\nval " << *obj_value << "\n";
  j->logMsg(ss.str());
  return TRUE;
}

Bool eval_grad_fj(Index n, Number* x, Bool new_x,
  Number* grad_f, UserDataPtr user_data)
{
  IpInterface * j = (IpInterface*)user_data;
  std::ostringstream ss;
  Eigen::VectorXd p = Eigen::VectorXd::Zero(n);
  ss << "grad: ";
  for (int ii = 0; ii<n; ii++){
    p[ii] = x[ii];
    ss << x[ii] << " ";
  }
  ss << "\n";
  if (j->p_cache.size() == p.size()){
    Eigen::VectorXd diff = j->p_cache - p;
    if (diff.squaredNorm()<1e-12){
      new_x = false;
    }
  }
  j->p_cache = p;

  //  std::cout<<"eval_grad_f new_x "<<" "<<has_grad_cache<<" "<<new_x<<"\n";
  if (!new_x && j->has_grad_cache){
    for (int ii = 0; ii<n; ii++){
      grad_f[ii] = j->grad_cache[ii];
      ss << grad_f[ii] << " ";
    }
    ss << "\n";
    j->logMsg(ss.str());
    return true;
  }

  j->logMsg(ss.str());
  ss.clear();
  j->objective->param = p;
  Eigen::VectorXd df = j->objective->df();
  for (Index i = 0; i < n; i++){
    grad_f[i] = df[i];
    ss << df[i] << " " ;
  }
  ss << "\n";
  j->logMsg(ss.str());

  j->has_grad_cache = true;
  j->grad_cache.resize(n);
  for (int ii = 0; ii<n; ii++){
    j->grad_cache[ii] = grad_f[ii];
  }
  if (new_x){
    j->has_f_cache = false;
  }

  return TRUE;
}

Bool eval_gj(Index n, Number* x, Bool new_x,
  Index m, Number* g, UserDataPtr user_data)
{
  return FALSE;
}

Bool eval_jac_gj(Index n, Number *x, Bool new_x,
  Index m, Index nele_jac,
  Index *iRow, Index *jCol, Number *values,
  UserDataPtr user_data)
{
  return FALSE;
}

Bool eval_hj(Index n, Number *x, Bool new_x, Number obj_factor,
  Index m, Number *lambda, Bool new_lambda,
  Index nele_hess, Index *iRow, Index *jCol,
  Number *values, UserDataPtr user_data)
{
  return FALSE;
}

Bool intermediate_cbj(Index alg_mod, Index iter_count, Number obj_value,
  Number inf_pr, Number inf_du, Number mu, Number d_norm,
  Number regularization_size, Number alpha_du,
  Number alpha_pr, Index ls_trials, UserDataPtr user_data)
{
  printf("callback in iteration %d\n", iter_count);
  IpInterface * j = (IpInterface*)user_data;
  j->objective->updateParam();
  return TRUE;
}

