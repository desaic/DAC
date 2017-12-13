#include "sparseMatrix.h"
#include "ARPACKSolver.h"
#include <iostream>
#include <fstream>

int example();

int main(int argc, char * argv[])
{
  if(argc<3){
    std::cout<<"eig stiffnessfile massfile\n";
    return 0;
  }
  SparseMatrix K(argv[1]);
  SparseMatrix M(argv[2]);

  printf("Computing linear modes using ARPACK: ...\n");

  int numModes = 10;
  double sigma = -1.0;
  int numLinearSolverThreads=1;
  double * frequenciesTemp = (double*) malloc (sizeof(double) * numModes);
  int DOFs = K.Getn();
  double * modesTemp = (double*) malloc
      (sizeof(double) * numModes * DOFs);

  ARPACKSolver generalizedEigenvalueProblem;
  int nconv = generalizedEigenvalueProblem.SolveGenEigShInv
      (&K, &M,
       numModes, frequenciesTemp,
       modesTemp, sigma, numLinearSolverThreads);

  std::cout<<"Arpack done with "<<nconv<<" modes\n";

  std::ofstream out("lambda.txt");
  out<<numModes<<"\n";
  for(int ii = 0; ii<numModes; ii++){
    out<<frequenciesTemp[ii]<<"\n";
  }
  out.close();

  out.open("modes.txt");
  for (int jj = 0; jj<DOFs; jj++){
    for(int ii = 0; ii<numModes; ii++){
      out<<modesTemp[ii*DOFs + jj]<<" ";
    }
    out<<"\n";
  }

  free(modesTemp);
  free(frequenciesTemp);

  return 0;
}
