#######################################################################################
# ukf.R
#
# Code to compute the unscented Kalman filter
#######################################################################################

#' Create the sigma nodes and weights for the unscented Kalman Filter
#' 
#' @param m The input mean
#' @param P The input covariance
#' @param n The dimension of the problem
#' @param alpha Scaling parameter
#' @param kappa Scaling parameter
#' @param betta Scaling parameter
#'   
#' @return A list containing: matirix of weights W in two rows, W^m and W^c; and
#'   the matrix of nodes X, where X(i) is in column i
#'   
sigma.weights <- function( m, P, n, alpha=1e-03, kappa=0, betta=2, quad=TRUE ){
  
  sqrt.P = t( chol( P ) )
  
  if( quad ){
    
    n.nodes <- 7
    
    nodes.1d <-  sqrt(2.0) * 
                  c( -2.651961356835233, -1.673551628767471, -0.8162878828589647, 0,
                      0.8162878828589647,  1.673551628767471,  2.651961356835233 )
    wts.1d <- c( 0.0009717812450995192, 0.05451558281912703, 0.4256072526101278, 
                 0.8102646175568073, 0.4256072526101278, 0.05451558281912703, 0.0009717812450995192 )
    wts.1d <- wts.1d / 1.77245385
    
    wts <- matrix( 0, 2, n.nodes ^ n )
    X <- matrix( 0, n, n.nodes ^ n )
    
    perms <- permutations( n.nodes, n, repeats.allowed = T )
    for( i in 1:(n.nodes^n) ){
      X[,i] <- m + sqrt.P * nodes.1d[ perms[ i, ] ]
      wts[,i] <- rep( prod( wts.1d[ perms[ i, ] ] ), 2 )
    }
    
  }else{
  
    lambda = alpha ^ 2 * ( n + kappa ) - n
        # The scaling parameter
    wts = .5 / ( n + lambda ) * matrix( 1, 2, 2*n+1 ) ;
        # Initialize the output matrix
    wts[ 1, 1 ] = lambda / ( n + lambda ) ; 
    wts[ 2, 1 ] = lambda / ( n + lambda ) + 1 - alpha ^ 2 + betta ;
        # Change the W_0 elements
    
    X <- matrix( 0, n, 2*n+1 )
        # Initialize the nodes
    X[,1] = m ;
    X[, 1 + 1:n ] = m + sqrt( n + lambda ) * sqrt.P ;
    X[, n + 1 + 1:n ] = m - sqrt( n + lambda ) * sqrt.P ;
        # The sigma points
    
  }
  
  out = list( wts=wts, X=X )
      # Create the output list
  return( out )
}


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
ukf.predict <- function( m, P, f, Q, n=length(m), alpha=1e-03, kappa=0, betta=2, quad=TRUE ){

  X.W <- sigma.weights( m, P, n, alpha, kappa, betta, quad )
  X <- X.W$X
  W <- X.W$wts
      # The integration nodes and weights (sigma points)
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
ukf.update <- function( m, P, g, R, y, n=length(m), alpha=1e-03, kappa=0, betta=2, quad=TRUE, n.mc=0 ){
  
  if( n.mc > 0 ){
    # Do Monte Carlo integration
    X <- m + mvrnorm( n.mc, m, P )
        # The sample of points
    Y <- matrix( apply( X, 1, g ), ncol=n.mc )
        # The matrix of measurements
    mu <- apply( Y, 1, mean )
    S <- var( t(Y) )
    C <- cov( X, matrix( t(Y), nrow=n.mc ))
  }else{
    X.W <- sigma.weights( m, P, n, alpha, kappa, betta, quad )
    X <- X.W$X
    W <- X.W$wts
        # The integration nodes and weights (sigma points)
    Y.hat <- matrix( apply( X, 2, g ), ncol=ncol(X) )
        # Evolve the sigma points through the observation equation
    mu <- Y.hat %*% W[1, ]
        # The predicted mean of the observation
    n.y <- length(mu)
        # The number of observations
    S <- R
    C <- matrix( 0, n, n.y )
        # Initialize the covariance matrices
    for( i in 1:ncol(W) ){
      S <- S + W[2,i] * ( Y.hat[,i] - mu ) %*% t( Y.hat[,i] - mu )
      C <- C + W[2,i] * ( X[,i] - m ) %*% t( Y.hat[,i] - mu )
    }    # Compute the covariances
  }
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
ukf.compute <- function( m0, P0, y, f, g, Q, R, n=length(m0), alpha=1e-03, kappa=0, betta=2, quad=TRUE, n.mc=0 ){

  n.y <- if( is.null( nrow(y) ) ) 1 else nrow(y)
  y <- matrix( y, nrow=n.y )
      # Matrix formatting
  Q <- matrix( Q, n, n )
  R <- matrix( R, n.y, n.y )
      # Make sure that y, Q and R are formatted as a matrix
  K <- ncol(y)
      # The number of points
  m <- matrix( 0, n, K )
  P <- array( 0, dim=c(dim(Q), K ) )  
#   m[,1] <- m0
#   P[,,1] <- P0
      # Initialize the mean and covariance matrices
  m.pred <- matrix( 0, n, K + 1 )
  P.pred <- array( 0, dim=c(dim(Q), K + 1 ) )
  m.pred[,1] <- m0
  P.pred[,,1] <- P0
      # The predictions for the mean and variance
  for( i in 1:K ){
    update <- ukf.update( m.pred[,i], matrix(P.pred[,,i], n, n ), g, R, y[,i], n, alpha, kappa, betta, quad, n.mc )
        # The updated mean and covariance
    m[,i] <- update$m
    P[,,i] <- update$P
        # Store
    pred <- ukf.predict( m[,i], matrix(P[,,i], n, n), f, Q, n, alpha, kappa, betta, quad )
        # The list of predicted values
    m.pred[,i+1] <- pred$m
    P.pred[,,i+1] <- pred$P
        # Store
  }
  return( list( m=m, P=P, m.pred=m.pred, P.pred=P.pred) )
}

rse <- function( ukf, x ){
# Computes the RSE of a ukf output
  rse <- sqrt( apply( ( ukf$m - x ) ^ 2, 2, sum ) )
  rse.pred <- sqrt( apply( ( matrix( ukf$m.pred[,-ncol(ukf$m.pred)], ncol=ncol(ukf$m.pred)-1) - x ) ^ 2, 2, sum ) )
  return( list( rse=rse, rse.pred=rse.pred ) )
}




