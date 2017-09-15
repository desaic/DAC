#include "Objective.hpp"
#include <iostream>

void Objective::getBounds(int idx, double &lb, double &ub)
{
  assert(idx<nx() && idx>0);
  lb = 1e-19;
  ub = 1e19;
}

void verifyObjective(double * x, Objective * obj)
{
  obj->setx(x);
  int nx = obj->nx();
  std::vector<double> df(nx, 0);
  int status = obj->df(df);
  std::vector<double> numdf(df.size(),0);
  double dx = 0.0001;
  double maxdiff = 0;

  for(int ii =0 ; ii<nx; ii++){
      double Ep,Em;

      x[ii] += dx;
      obj->setx(x);
      Ep = obj->f(status);

      x[ii] -= 2*dx;
      obj->setx(x);
      Em = obj->f(status);

      x[ii] += dx;
      obj->setx(x);

      double diff = 0.5*(Ep-Em)/dx;
      numdf[ii] = diff;
      maxdiff = std::max(maxdiff, (double)std::abs(diff - df[ii]));
      std::cout<<ii<<" "<<df[ii]<<" "<<numdf[ii]<<" "<<diff<<"\n";
  }
  std::cout<<"force maxdiff "<<maxdiff<<"\n";

  maxdiff = 0;
  int maxIdx =0 ;
  //check hessian
  Eigen::SparseMatrix<double> ddf = obj->d2f();
  for(int ii = 0; ii<nx; ii++){

    x[ii] += dx;
    obj->setx(x);
    std::vector<double> dfp(nx,0);
    obj->df(dfp);

    x[ii] -= 2*dx;
    obj->setx(x);

    std::vector<double> dfm(nx,0);
    obj->df(dfm);

    for (Eigen::SparseMatrix<double>::InnerIterator it(ddf, ii); it; ++it){
      int jj = it.row();
      dfp[jj] -= dfm[jj];
      dfp[jj] /= 2*dx;
      std::cout<<ii<<" "<<jj<<" "<<it.value()<<" "<<dfp[jj]<<"\n";
      float diff = (float)std::abs(dfp[jj]-it.value());
      if(diff>maxdiff){
        maxdiff = diff;
        maxIdx = ii;
      }
    }

    x[ii] += dx;
    obj->setx(x);
  }
  std::cout<<"stiffness maxdiff "<<maxdiff<<" "<<maxIdx<<"\n";

} 
