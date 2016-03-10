/***********************************************************************************
 * ukf.hpp
 * 
 * Interface to ukf.cpp
 * 
 * 09mar2016
 * Philip Barrett, Chicago
 * 
 ***********************************************************************************/

#ifndef UKF_HPP
#define UKF_HPP

#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

arma::mat sigma_pts( arma::vec m, arma::mat P, int n, double alpha, double kappa, double betta ) ;
    // Computes the sigma points used in integration

#endif