#include "ArrayUtil.hpp"
#include "NewtonMin.hpp"
#include "Objective.hpp"
#include "Timer.hpp"
#include <Eigen/SparseCholesky>
#include <iostream>
#include <fstream>

int NewtonMin::minimize(std::vector<double> & x)
{
  //verifyObjective(x.data(), obj);

  std::vector<double> df(x.size(),0);
  std::vector<double> x0;
  std::vector<double> delta(x.size(),0);
  Eigen::SparseMatrix<double> d2f;
  int status;
  int ret = 0;
  Timer timer;
  double n0 = 0, n1 = 0, n2 = 0;
  int input;
  for(int iter = 0; iter<maxIters; iter++){
    //std::cout << "Newton min iter " << iter << "\n";
    float steplen=1.0f;
    obj->setx(x.data());
    //timer.startWall();
    //std::cout << "assemble stiffness \n";
    //std::cin >> input;
    d2f = obj->d2f();
//#ifdef _DEBUG
    //Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > DEC(d2f);
    //if (DEC.info() == Eigen::NumericalIssue){
    //  throw std::runtime_error("BAD!");
    //}
//#endif
    //timer.endWall();
    //std::cout<<"stiffness total time "<<timer.getSecondsWall()<<"\n";
    //std::cout << "lin solve\n";
    //std::cin >> input;
    status = obj->df(df);
    n1 = infNorm(df);
    if(iter==0){
      n0 = n1;
    }
    if(iter>0 &&( (n1/n0) < df_rel_tol || n1 < df_abs_tol ) ){
      return ret;
    }
    timer.startWall();
    double residual;
    int iters=0;
    
    solver->solve(d2f, df, delta, residual, iters);
    //std::cout << "lin done\n";
    //std::cin >> input;

    timer.endWall();
    std::cout<<"lin solve time "<<timer.getSecondsWall()<<"\n";
    //std::cout<<"NewtonMin lin residual "<<residual<<"\n";
    if(residual>residual_tol){
      //    write_matlab(std::cout, "K", d2f);
//      std::cout<<"NewtonMin large lin residual "<<residual<<"\n";
//      ret = -1;
    }
    std::cout<<"newton "<< iter<<" "<<n1<<"\n";
    x0 = x;
    double E0 = obj->f(status);
    while(1){
      x = x0;
      //minimization steps in the negative gradient direction
      add_scaled(x, -steplen, delta);
      obj->setx(x.data());
      status = obj->df(df);
      int statusE = 0;
      double E1 = obj->f(statusE);
      n2 = infNorm(df);
      std::cout << "E1 E0 " << E1 << " " << E0 << "\n";
      std::cout<<"h n1 n2 " <<steplen<<" "<<n1<<" "<<n2<<"\n";
      if( (E1>E0 && n2>n1) || status<0){
        steplen = 0.5f * steplen;
        if(steplen<1e-6){
          std::cout << "NewtonMin.cpp cannot find step size that reduces gradient norm\n";
          ret = -2;
          break;
        }
      }else{
         break;
      }
    }
    if(ret<0){
      return ret;
    }
    n1 = n2;

    //if(iter==0){
    //  std::cout<<iter<<" norm "<<n0<<"\n";
    //}
    //std::cout << (iter+1) << " norm " << n1<<"\n";

  }
  //std::cout << "=========\n";
  if(  (n1/n0) < df_rel_tol || n1 < df_abs_tol  ){
    //std::cout<<"h n1 n2 " <<n1<<" "<<n0<<"\n";
    return ret;
  }
//  std::cout<<"Newton min error "<<ret<<"\n";
  //did not converge within given number of steps.
  return -3;
}
