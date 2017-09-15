#ifndef OBJECTIVE_HPP
#define OBJECTIVE_HPP

#include <vector>
#include <Eigen/Sparse>

class Objective{
public:
  virtual double f(int & staus)=0;
  virtual int df(std::vector<double> & grad)=0;
  virtual Eigen::SparseMatrix<double> d2f()=0;
  virtual void setx(double * x)=0;

  ///@brief number of variables in the state.
  virtual int nx()=0;

  virtual void getBounds(int idx, double & lb, double &ub);

  virtual ~Objective(){}
}; 

////==============================================
////numerical differencing test
void verifyObjective(double *x, Objective *obj);

#endif
