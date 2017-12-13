/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.0                               *
 *                                                                       *
 * "sparseMatrix" library , Copyright (C) 2007 CMU, 2009 MIT, 2013 USC   *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Jernej Barbic                                            *
 * http://www.jernejbarbic.com/code                                      *
 *                                                                       *
 * Research: Jernej Barbic, Fun Shing Sin, Daniel Schroeder,             *
 *           Doug L. James, Jovan Popovic                                *
 *                                                                       *
 * Funding: National Science Foundation, Link Foundation,                *
 *          Singapore-MIT GAMBIT Game Lab,                               *
 *          Zumberge Research and Innovation Fund at USC                 *
 *                                                                       *
 * This library is free software; you can redistribute it and/or         *
 * modify it under the terms of the BSD-style license that is            *
 * included with this library in the file LICENSE.txt                    *
 *                                                                       *
 * This library is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

/*
  Example for how to use the SparseMatrix class and the conjugate gradient solver.
  The CG solver is in the "sparseSolvers" folder.
  See also sparseMatrix.h and CGSolver.h.
*/

#include "sparseMatrix.h"
#include "CGSolver.h"
#include "ARPACKSolver.h"
#include <iostream>
#include <fstream>

int example()
{

  int n = 100; // matrix dimension

  SparseMatrixOutline * outline = new SparseMatrixOutline(n); // n x n sparse matrix
  for(int row=0; row<n; row++){
    // here we are also setting the values of non-zero entries (or we could also do that later with the SparseMatrix object)
    if (row - 1 >= 0)
      outline->AddEntry(row, row-1, 1.0);
    outline->AddEntry(row, row, 3.0);
    if (row + 1 < n)
      outline->AddEntry(row, row+1, 1.0);
  }

  outline->Save("A.txt");
  SparseMatrix A(outline);
  delete(outline);

  outline = new SparseMatrixOutline(n); // n x n sparse matrix
  for(int row=0; row<n; row++){
    outline->AddEntry(row, row, 1.0);
  }

  SparseMatrix M(outline);
  delete(outline);

  printf("Computing linear modes using ARPACK: ...\n");

  int numDesiredModes = 10;
  double sigma = 0;
  int numLinearSolverThreads=3;
  double * frequenciesTemp = (double*) malloc (sizeof(double) * numDesiredModes);
  int numRetainedDOFs = A.Getn();
  double * modesTemp = (double*) malloc
      (sizeof(double) * numDesiredModes * numRetainedDOFs);

  ARPACKSolver generalizedEigenvalueProblem;
  int nconv = generalizedEigenvalueProblem.SolveGenEigShInv
      (&A, &M,
       numDesiredModes, frequenciesTemp,
       modesTemp, sigma, numLinearSolverThreads);

  std::cout<<"Arpack done with "<<nconv<<" modes\n";

  for(int ii = 0; ii<numDesiredModes; ii++){
    for (int jj = 0; jj<numRetainedDOFs; jj++){
      std::cout<<modesTemp[ii*numRetainedDOFs + jj]<<" ";
    }
    std::cout<<"\n";
  }

  free(modesTemp);
  free(frequenciesTemp);

  return 0;
}
