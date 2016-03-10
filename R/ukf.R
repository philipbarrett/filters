#######################################################################################
# ukf.R
#
# Code to compute the unscented Kalman filter
#######################################################################################

#' Predict the state for the unscented Kalman Filter
#' 
#' @param m The previous mean
#' @param P The previous covariance
#' @param f The law of motion for the state
#' @param Q The state covariance
#' @param n The dimension of the problem
#' 
#' @return The mean and covariance of the predicted state
#' 
ukf.predict <- function( m, P, f, Q, n=length(m), alpha=1e-03, kappa=0, betta=2 ){

  X <- sigma_pts( m, P, n, alpha, kappa, betta )
      # The integration nodes (sigma points)
  W <- sigma_weights( n, alpha, kappa, betta )
      # The integration weights
  X.hat <- matrix( apply( X, 2, f ), nrow=n )
      # Evolve the sigma points through the state law of motion
  m.new <- X.hat %*% W[1, ]
      # The predicted mean
  P.new <- Q
  for( i in 1:(2*n+1))
    P.new <- P.new + W[2,i] * ( X.hat[,i] - m.new ) %*% t( X.hat[,i] - m.new )
      # The predicted covariance
  return( list( m=m.new, P=P.new ) )
}

#' Predict the state for the unscented Kalman Filter
#' 
#' @param m The predicted mean
#' @param P The predicted covariance
#' @param g The measurement equation
#' @param R The measurement covariance
#' @param n The dimension of the problem
#' 
#' @return The mean and covariance of the predicted state
#' 
ukf.update <- function( m, P, g, R, y, n=length(m), alpha=1e-03, kappa=0, betta=2 ){
  
  X <- sigma_pts( m, P, n, alpha, kappa, betta )
      # The integration nodes (sigma points)
  W <- sigma_weights( n, alpha, kappa, betta )
      # The integration weights
  Y.hat <- matrix( apply( X, 2, g ), ncol=2*n+1 )
      # Evolve the sigma points through the observation equation
  mu <- Y.hat %*% W[1, ]
      # The predicted mean of the observation
  n.y <- length(mu)
      # The number of observations
  S <- R
  C <- matrix( 0, n, n.y )
      # Initialize the covariance matrices
  for( i in 1:(2*n+1) ){
    S <- S + W[2,i] * ( Y.hat[,i] - mu ) %*% ( Y.hat[,i] - mu )
    C <- C + W[2,i] * ( X[,i] - m ) %*% t( Y.hat[,i] - mu )
  }   # Compute the covariances
  K <- C %*% solve( S )
      # The gain
  m.new <- m + K %*% ( y - mu )
      # The updated mean
  P.new <- P - K %*% S %*% t( K )
      # The updated covariance
  return( list( m=m.new, P=P.new ) )
}

#' Compute the unscented Kalman Filter
#' 
#' @param m0 The initial mean
#' @param P0 The initial covariance
#' @param y The sequence of signals
#' @param f The law of motion for the state
#' @param g The measurement equation
#' @param Q The state shock covariance
#' @param R The measurement covariance
#' @param n The dimension of the problem
#' 
#' @return The mean and covariance of the predicted state at each point in time
#' 
ukf.compute <- function( m0, P0, y, f, g, Q, R, n=length(m), alpha=1e-03, kappa=0, betta=2 ){

  if( is.null( nrow(y) ) ) y <- matrix( y, nrow=1 )
  Q <- matrix( Q )
  R <- matrix( R )
      # Make sure that y, Q and R are formatted as a matrix
  K <- ncol(y)
      # The number of points
  m <- matrix( 0, n, K + 1 )
  P <- array( 0, dim=c(dim(Q), K+1 ) )  
  m[,1] <- m0
  P[,,1] <- P0
      # Initialize the mean and covariance matrices
  m.pred <- matrix( 0, n, K )
  P.pred <- array( 0, dim=c(dim(Q), K ) )
      # The predictions for the mean and variance
  for( i in 1:K ){
    pred <- ukf.predict( m[,i], matrix(P[,,i]), f, Q, n, alpha, kappa, betta )
        # The list of predicted values
    m.pred[,i] <- pred$m
    P.pred[,,i] <- pred$P
        # Store
    update <- ukf.update( pred$m, pred$P, g, R, y[,i], n, alpha, kappa, betta )
        # The updated mean and covariance
    m[,i+1] <- update$m
    P[,,i+1] <- update$P
        # Store
  }
  return( list( m=m, P=P, m.pred=m.pred, P.pred=P.pred) )
}