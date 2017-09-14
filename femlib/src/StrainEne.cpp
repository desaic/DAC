#include "StrainEne.hpp"
#include <iostream>

StrainEne::~StrainEne(){}

#include "StrainCorotLin.hpp"
#include "StrainLin.hpp"
#include "StrainEneNeo.hpp"
#include "StrainEneANeo.hpp"

///@param E array of 2 doubles. E and nu
///@param mu array of size 2. mu and lambda.
void Young2Lame(double * E, double * mu){
  mu[0] = E[0] / (2*(1+E[1]));
  mu[1] = E[0] * E[1]/((1+E[1])*(1-2*E[1]));
}

std::vector < StrainEne* > loadMaterials(std::istream & in, std::vector<double> & densities)
{
  int nMat = 0;
  in>>nMat;
  std::vector<StrainEne*> materials(nMat, 0);
  for(int ii = 0; ii<nMat; ii++){
    int matType=0;
    in >>matType;
    int nParam = 0;
    in>>nParam;
    std::vector<double> param(nParam, 0);
    for(int jj = 0; jj<nParam; jj++){
      in>>param[jj];
    }
    StrainEne * ene = 0;
    float density = 0;
    bool knownMat = true;

    switch(matType){
    case StrainEne::LIN:
      ene = new StrainLin();
      ene->param.resize(2);
      Young2Lame( param.data(), ene->param.data() );
      break;

    case StrainEne::COROT:
      ene = new StrainCorotLin();
      ene->param.resize(2);
      Young2Lame( param.data(), ene->param.data() );
      break;

    case StrainEne::NEO:
      ene = new StrainEneNeo();
      ene->param.resize(2);
      Young2Lame( param.data(), ene->param.data() );
      break;
    case StrainEne::ANEO:
  	  ene = new StrainEneANeo();
      ene->param=param;
      Young2Lame( param.data(), ene->param.data() );
      break;
    default:
      std::cout<<"Unrecognized material type "<<matType<<"\n";
      knownMat = false;
      break;
    }
    if (knownMat){
      in >> density;
    }
    densities.push_back(density);
    materials[ii] = ene;
  }
  return materials;
}
