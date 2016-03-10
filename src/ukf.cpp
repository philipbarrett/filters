/***********************************************************************************
 * ukf.hpp
 * 
 * Code to compute the unscented Kalman filter
 * 
 * 09mar2016
 * Philip Barrett, Chicago
 * 
 ***********************************************************************************/

#include "ukf.hpp"

// [[Rcpp::export]]
arma::mat sigma_weights( int n, double alpha=1e-03, double kappa=0, double betta=2 ){
// Calculates the weights applied to the sigma points in two rows; W^m and then W^c
  
  double lambda = pow( alpha, 2 ) * ( n + kappa ) - n ;
      // The scaling parameter
  mat out = .5 / ( n + lambda ) * ones( 2, 2*n+1 ) ;
      // Initialize the output matrix
  out( 0, 0 ) = lambda / ( n + lambda ) ; 
  out( 1, 0 ) = lambda / ( n + lambda ) + 1 - pow(alpha, 2) + betta ;
      // Change the W_0 elements
  return out ;
}

// [[Rcpp::export]]
arma::mat sigma_pts( arma::vec m, arma::mat P, int n, double alpha=1e-03, double kappa=0, 
                      double betta=2 ){
// Computes the sigma points used in integration, each column is one of the sigma points
  double lambda = pow( alpha, 2 ) * ( n + kappa ) - n ;
      // The scaling parameter
  mat out = zeros( n, 2*n+1 ) ;
      // Initialize output
  mat sqrt_P = chol( P, "lower" ) ;
      // The Cholesky decomposition: sqrt_P * sqrt_P.t() = P
  mat sqrt_P_t = sqrt_P.t() ;
      // Transpose (means that the X^i can be formed as columns more easily)
  out.col(0) = m ;
  out.cols( 1, n ) = m * ones<rowvec>(n) + sqrt( n + lambda ) * sqrt_P_t ;
  out.cols( n+1, 2*n ) = m * ones<rowvec>(n) - sqrt( n + lambda ) * sqrt_P_t ;
      // The sigma points
  return out ;
}

// [[Rcpp::export]]
arma::mat chol_test( arma::mat M ){
// Testing the Cholseky decomp
  mat out = chol( M, "lower" ) ;
  return out ;
}