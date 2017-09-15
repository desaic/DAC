#ifndef RUN_LBFGS_HPP
#define RUN_LBFGS_HPP

#include "JumperObjective.hpp"
#include "lbfgs.h"

lbfgsfloatval_t _evaluate(
  void *instance,
  const lbfgsfloatval_t *x,
  lbfgsfloatval_t *g,
  const int n,
  const lbfgsfloatval_t step
  );

int _progress(
  void *instance,
  const lbfgsfloatval_t *x,
  const lbfgsfloatval_t *g,
  const lbfgsfloatval_t fx,
  const lbfgsfloatval_t xnorm,
  const lbfgsfloatval_t gnorm,
  const lbfgsfloatval_t step,
  int n,
  int k,
  int ls);

int runLBFGS(const ConfigFile & config,
  JumperObjective &objective);

#endif