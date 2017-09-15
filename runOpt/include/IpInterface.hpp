// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id: MyNLP.hpp 1861 2010-12-21 21:34:47Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-11-05

#ifndef IPINTERFACE_HPP
#define IPINTERFACE_HPP

#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include "IpStdCInterface.h"

class RealFun;

class IpInterface
{
public:
  /** default constructor */
  IpInterface();

  /** default destructor */
  virtual ~IpInterface();
  int init();
  int run();
  void logMsg(const std::string & msg)
  {
    if (logger != 0 && logger->good()){
      std::cout << msg << "\n";
      (*logger) << msg;
      logger->flush();
    }
  }

  RealFun * objective;
  std::vector<Number> x_L, x_U;
 
  float f_cache;
  std::vector<double> grad_cache;
  bool has_f_cache;
  bool has_grad_cache;
  Eigen::VectorXd p_cache;

  std::ostream * logger;

private:
  IpoptProblem nlp;
};


#endif
