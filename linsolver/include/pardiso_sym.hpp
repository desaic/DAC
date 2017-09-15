#ifndef PARDISO_SYM_HPP
#define PARDISO_SYM_HPP

struct PardisoState{
  /* Internal solver memory pointer pt,                  */
  /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
  /* or void *pt[64] should be OK on both architectures  */
  void    *pt[64];

  /* Pardiso control parameters. */
  int      iparm[64];
  double   dparm[64];
  int      maxfct, mnum, msglvl;
  /* Auxiliary variables. */
  int      mtype;        /* Real symmetric matrix */
};

void pardisoInit(PardisoState * state);
int pardisoSymbolicFactorize(int * ia, int * ja, int n, PardisoState *state);
int pardisoNumericalFactorize(int * ia, int * ja, double* val, int n, PardisoState *state);
int pardisoBackSubstitute(int * ia, int * ja, double* val, int n, double* x, double * b, PardisoState *state);
int pardisoSolve(int * ia, int * ja, double* val, int n, double* x, double * b, PardisoState *state);
void pardisoFree(int * ia, int * ja, int n, PardisoState *state);

#endif 
