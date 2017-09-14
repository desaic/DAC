#include "StepperGrad.hpp"
#include "ElementMesh.hpp"
#include "femError.hpp"
#include "ArrayUtil.hpp"

#include <Eigen/Dense>
#include <vector>
#include <iostream>

StepperGrad::StepperGrad():h(0.01f),force_L2tol(1e-4f) {}

int StepperGrad::oneStep()
{
  int dim = m->dim;
  std::vector<Eigen::Vector3d> force = m->getForce();
  for(unsigned int ii = 0;ii<force.size();ii++){
    for(int jj = 0; jj<dim; jj++){
    if(m->fixedDof[dim * ii + jj]){
      force[ii][jj] = 0;
    }
    }
  }
  double E = m->getEnergy();
  double totalMag = 0;
  for(unsigned int ii = 0;ii<force.size();ii++){
    totalMag += force[ii].transpose()*force[ii];
  }
  if(totalMag<force_L2tol){
    return 0;
  }

  std::vector<Eigen::Vector3d> x0 = m->x;
  while(1){
    m->x=x0;
    addmul(m->x, h, force);
    fem_error = 0;
    double E1 = m->getEnergy();
    std::cout << "h " << h << "\n";
    std::cout << "E " << E1 << " "<< E<<"\n";
    if(E1>E || fem_error){
      fem_error = 0;
      h = 0.5f* h;
      std::cout<<"h "<<h<<"\n";
    }else{
      h=1.1f*h;
      break;
    }
  }
  return 0;
}
