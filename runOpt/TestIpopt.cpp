#include "TestIpopt.hpp"
#include "IpStdCInterface.h"
#include "IpInterface.hpp"
#include "RealFun.hpp"
#include <fstream>
/* Ipopt Function Declarations */
Bool eval_f(Index n, Number* x, Bool new_x,
  Number* obj_value, UserDataPtr user_data);

Bool eval_grad_f(Index n, Number* x, Bool new_x,
  Number* grad_f, UserDataPtr user_data);

Bool eval_g(Index n, Number* x, Bool new_x,
  Index m, Number* g, UserDataPtr user_data);

Bool eval_jac_g(Index n, Number *x, Bool new_x,
  Index m, Index nele_jac,
  Index *iRow, Index *jCol, Number *values,
  UserDataPtr user_data);

Bool eval_h(Index n, Number *x, Bool new_x, Number obj_factor,
  Index m, Number *lambda, Bool new_lambda,
  Index nele_hess, Index *iRow, Index *jCol,
  Number *values, UserDataPtr user_data);

Bool intermediate_cb(Index alg_mod, Index iter_count, Number obj_value,
  Number inf_pr, Number inf_du, Number mu, Number d_norm,
  Number regularization_size, Number alpha_du,
  Number alpha_pr, Index ls_trials, UserDataPtr user_data);

/* This is an example how user_data can be used. */
struct MyUserData
{
  Number g_offset[2]; /* This is an offset for the constraints.  */
};

class TestObjective:public RealFun{
public:
  double lambda[2];
  double rho;
  double g[2];
  TestObjective() :rho(1e4){
    lambda[0] = 0;
    lambda[1] = 0;
    g[0] = 0; 
    g[1] = 0;
  }
  void updateParam(){
    f();
    for (int i = 0; i<2; i++){
      lambda[i] = rho * std::max(g[i] + lambda[i] / rho, 0.0);
    }
  }

  double f(){
    Eigen::Vector2d v0;
    double w0 = 51;
    double clearWall = 1e-2;
    v0 << 0.2, 1.6;//0.23 1.8
    v0[0] += param[0] + 15 * param[1];
    v0[1] += -10 * param[0] + 100 * param[1];
    w0 += -param[0] + 100 * param[1];

    Eigen::Vector2d grav;
    grav << 0, -9.8;
    
    double dt = 1e-4;
    int nStep = (int)(std::abs(2 * v0[1] / grav[1]/dt) + 1);
    std::vector<Eigen::Vector2d> traj;
    Eigen::Vector2d wall;
    wall << 0.05, 0.1;
    bool hitwall = false;
    Eigen::Vector2d hitPt;
    double eps = 1e-6;
    double T = 0;
    for (int i = 0; i < nStep; i++){
      double t = i*dt;
      Eigen::Vector2d x = v0 * t + 0.5 * t * t * grav;
      if (traj.size() > 0){
        Eigen::Vector2d x0 = traj[traj.size() - 1];
        if (x0[0]<=wall[0] && x[0]>=wall[0]){
          if ((x0[1] < wall[1] || x[1] < wall[1])){
            hitwall = true;
            hitPt = x0;
            break;
          }
        }
      }
      if (x[1] < -eps){
        break;
      }
      T = t;
      traj.push_back(x);
    }
    double theta = T * w0;
    double pi = 3.14159265359;
    int nflips = theta/(2*pi);
    double phi = theta - nflips * 2 * pi;
    double val = 0.5 * phi * phi;
    Eigen::Vector2d lastx = traj[traj.size() - 1];
    if (hitwall){
      g[0] = wall[1] - (lastx[1] - clearWall);
    }
    else{
      g[0] = 0;
    }
    for (int i = 0; i < 2; i++){
      double ghat = std::max(g[i] + lambda[i] / rho, 0.0);
      val += 0.5*rho * ghat * ghat - 0.5 * lambda[i] * lambda[i] / rho;
    }
    return val;
  }
  
  ///@brief compute gradient. Not implemented by default.
  ///The return value should have the same dimensions as the parameters.
  virtual Eigen::VectorXd df(){
    Eigen::VectorXd grad = Eigen::VectorXd::Zero(param.size());
    //std::cout << "grad\n";
    grad[0] = 2 * (param[0]-0.5);
    grad[1] = 4 * (param[1] -0.8);
    double h = 2e-5;
    for (int i = 0; i < grad.size(); i++){
      double vp, vm;
      param[i] += h;
      vp = f();
      param[i] -= 2 * h;
      vm = f();
      param[i] += h;
      grad[i] = (vp - vm) / (2 * h);
    }
    return grad;
  }
};

void testIpInterface()
{
  IpInterface ip;
  TestObjective * obj = new TestObjective();
  ip.logger = new std::ofstream("debugopt.txt");
  Eigen::VectorXd x0(2);
  x0[0] = 0;
  x0[1] = 0;
  obj->setParam(x0);
  ip.x_L=std::vector<double>(2,-0.005);
  ip.x_U = std::vector<double>(2, 0.005);
  ip.objective = obj;
  ip.init();
  ip.run();
}

void testIpopt()
{
  Index n = -1;                          /* number of variables */
  Index m = -1;                          /* number of constraints */
  std::vector<Number> x_L,x_U;         /* lower bounds on x upper bounds on x */
  IpoptProblem nlp = NULL;             /* IpoptProblem */
  enum ApplicationReturnStatus status; /* Solve return code */
                 /* starting point and solution vector */
  std::vector<Number>x, mult_x_L, mult_x_U;/* lower bound multipliers
                                       at the solution upper bound multipliers
                                       at the solution */
  Number obj;                          /* objective value */
  Index i;                             /* generic counter */

  /* Number of nonzeros in the Jacobian of the constraints */
  Index nele_jac = 0;
  /* Number of nonzeros in the Hessian of the Lagrangian (lower or
  upper triangual part only) */
  Index nele_hess = 0;
  /* indexing style for matrices */
  Index index_style = 0; /* C-style; start counting of rows and column
                         indices at 0 */

  /* our user data for the function evalutions. */
  struct MyUserData user_data;

  /* set the number of variables and allocate space for the bounds */
  n = 2;
  x_L.resize(n);
  x_U.resize(n);
  /* set the values for the variable bounds */
  for (i = 0; i<n; i++) {
    x_L[i] = 2.0;
    x_U[i] = 3.0;
  }

  /* set the number of constraints and allocate space for the bounds */
  m = 0;

  /* create the IpoptProblem */
  nlp = CreateIpoptProblem(n, x_L.data(), x_U.data(), m, NULL, NULL, nele_jac, nele_hess,
    index_style, &eval_f, &eval_g, &eval_grad_f,
    &eval_jac_g, &eval_h);

  /* Set some options.  Note the following ones are only examples,
  they might not be suitable for your problem. */
  AddIpoptNumOption(nlp, "tol", 1e-7);
  AddIpoptStrOption(nlp, "mu_strategy", "adaptive");
  AddIpoptStrOption(nlp, "hessian_approximation", "limited-memory");
  //AddIpoptStrOption(nlp, "output_file", "ipopt.out");

  /* allocate space for the initial point and set the values */
  x.resize(n);
  x[0] = 2.5;
  x[1] = 2.5;

  mult_x_L.resize(n);
  mult_x_U.resize(n);

  /* Initialize the user data */
  user_data.g_offset[0] = 0.;
  user_data.g_offset[1] = 0.;

  /* Set the callback method for intermediate user-control.  This is
  * not required, just gives you some intermediate control in case
  * you need it. */
  SetIntermediateCallback(nlp, intermediate_cb);

  /* solve the problem */
  status = IpoptSolve(nlp, x.data(), NULL, &obj, NULL, NULL, NULL, &user_data);

  if (status == Solve_Succeeded) {
    printf("\n\nSolution of the primal variables, x\n");
    for (i = 0; i<n; i++) {
      printf("x[%d] = %e\n", i, x[i]);
    }

    printf("\n\nSolution of the bound multipliers, z_L and z_U\n");
    for (i = 0; i<n; i++) {
      printf("z_L[%d] = %e\n", i, mult_x_L[i]);
    }
    for (i = 0; i<n; i++) {
      printf("z_U[%d] = %e\n", i, mult_x_U[i]);
    }

    printf("\n\nObjective value\n");
    printf("f(x*) = %e\n", obj);
  }
  else {
    printf("\n\nERROR OCCURRED DURING IPOPT OPTIMIZATION.\n");
  }

  /* free allocated memory */
  FreeIpoptProblem(nlp);
}

/* Function Implementations */
Bool eval_f(Index n, Number* x, Bool new_x,
  Number* obj_value, UserDataPtr user_data)
{
  *obj_value = 5 * x[0] * x[0] + 100 * x[1] * x[1] + 5;
  std::cout << "val " << x[0] << " " << x[1] << "\n" << *obj_value<<"\n";
  return TRUE;
}

Bool eval_grad_f(Index n, Number* x, Bool new_x,
  Number* grad_f, UserDataPtr user_data)
{
  
  grad_f[0] = 10 * x[0];
  grad_f[1] = 200 * x[1];
  std::cout << "grad " << x[0] << " " << x[1] << "\n" << grad_f[0]<<" "<<grad_f[1] << "\n";
  return TRUE;
}

Bool eval_g(Index n, Number* x, Bool new_x,
  Index m, Number* g, UserDataPtr user_data)
{
  struct MyUserData* my_data = (MyUserData*)user_data;
   
  return FALSE;
}

Bool eval_jac_g(Index n, Number *x, Bool new_x,
  Index m, Index nele_jac,
  Index *iRow, Index *jCol, Number *values,
  UserDataPtr user_data)
{
  
  return FALSE;
}

Bool eval_h(Index n, Number *x, Bool new_x, Number obj_factor,
  Index m, Number *lambda, Bool new_lambda,
  Index nele_hess, Index *iRow, Index *jCol,
  Number *values, UserDataPtr user_data)
{
  return FALSE;
}

Bool intermediate_cb(Index alg_mod, Index iter_count, Number obj_value,
  Number inf_pr, Number inf_du, Number mu, Number d_norm,
  Number regularization_size, Number alpha_du,
  Number alpha_pr, Index ls_trials, UserDataPtr user_data)
{
  printf("Testing intermediate callback in iteration %d\n", iter_count);

  return TRUE;
}

