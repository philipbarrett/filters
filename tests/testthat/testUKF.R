library(filters)
library(FKF)
library(MASS)

context("Comparing unscented to standard Kalman Filter")

sim.data <- function( x.0, f, g, R, Q, K ){
# Creates data to test the unscented Kalman Filter
  n.x <- nrow(Q)
  n.y <- nrow(R)
      # Dimensions of the state and measurement
  x <- matrix(0, n.x, K)
  y <- matrix(0, n.y, K)
      # Initialize the outputs
  x[,1] <- x.0
  y[,1] <- g( x[,1] ) + mvrnorm( 1, rep(0, nrow(R)), R )
      # The first columns
  for( i in 2:K ){
    x[,i] <- f( x[,i-1] ) + mvrnorm( 1, rep(0, nrow(Q)), Q )
    y[,i] <- g( x[,i] ) + mvrnorm( 1, rep(0, nrow(R)), R )
  }   # Create the state and measurment equations
  return( list( x=x, y=y ) )
}

test_that("A one-dimensional linear example", {

  K <- 100
      # Number of points
  FF <- matrix( .9, 1, 1 )
      # An AR(1)
  GG <- matrix( 1, 1, 1 )
      # Unbiased measurement
  f <- function( x ) FF %*% x
  g <- function( x ) GG %*% x
      # The evolutions
  Q <- matrix( 1, 1 , 1)   # State
  R <- matrix( 1, 1 , 1)   # Obs
      # Unit variances
  sim <- sim.data( 0, f, g, R, Q, K )
      # Generate the data
  m.0 <- 0
  P.0 <- matrix( 10, 1, 1 )
      # Initialize the mean and covariance
  base <- fkf( m.0, P.0, 0, 0, FF, GG, Q, R, sim$y )
      # From the FKF package
  alt <- ukf.compute( m.0, P.0, sim$y, f, g, Q, R, quad=F )
      # the unscented filter
  expect_equal( base$att, alt$m )
  expect_equal( base$at, alt$m.pred )
  expect_equal( base$Ptt, alt$P )
  expect_equal( base$Pt, alt$P.pred )
  
  ## Make observation noise bigger
  R <- matrix( 2, 1 , 1)
  sim <- sim.data( 0, f, g, R, Q, K )
  base <- fkf( m.0, P.0, 0, 0, FF, GG, Q, R, sim$y )
      # From the FKF package
  alt <- ukf.compute( m.0, P.0, sim$y, f, g, Q, R, quad=F )
      # the unscented filter
  expect_equal( base$att, alt$m )
  expect_equal( base$at, alt$m.pred )
  expect_equal( base$Ptt, alt$P )
  expect_equal( base$Pt, alt$P.pred )
  
  ## Make state noise bigger
  Q <- matrix( 4, 1 , 1)
  sim <- sim.data( 0, f, g, R, Q, K )
  base <- fkf( m.0, P.0, 0, 0, FF, GG, Q, R, sim$y )
      # From the FKF package
  alt <- ukf.compute( m.0, P.0, sim$y, f, g, Q, R, quad=F )
      # the unscented filter
  expect_equal( base$att, alt$m )
  expect_equal( base$at, alt$m.pred )
  expect_equal( base$Ptt, alt$P )
  expect_equal( base$Pt, alt$P.pred )
  
  ## Change observation equation 
  GG <- matrix( 2, 1, 1 )
  hh <- c( .5 )
      # Biased measurement
  g <- function( x ) hh + GG %*% x
      # The new measurement equation
  sim <- sim.data( 0, f, g, R, Q, K )
  base <- fkf( m.0, P.0, 0, hh, FF, GG, Q, R, sim$y )
      # From the FKF package
  alt <- ukf.compute( m.0, P.0, sim$y, f, g, Q, R, quad=F )
      # the unscented filter
  expect_equal( base$att, alt$m )
  expect_equal( base$at, alt$m.pred )
  expect_equal( base$Ptt, alt$P )
  expect_equal( base$Pt, alt$P.pred )
  
  ## Unit root in state
  FF <- matrix( 1, 1, 1 )
      # Unit root
  sim <- sim.data( 0, f, g, R, Q, K )
  base <- fkf( m.0, P.0, 0, hh, FF, GG, Q, R, sim$y )
      # From the FKF package
  alt <- ukf.compute( m.0, P.0, sim$y, f, g, Q, R, quad=F )
      # the unscented filter
  expect_equal( base$att, alt$m )
  expect_equal( base$at, alt$m.pred )
  expect_equal( base$Ptt, alt$P )
  expect_equal( base$Pt, alt$P.pred )
  
})

test_that("A two-dimensional linear example", {
# 2D tracking example taken from Sarkka Ch.4
  
  K <- 100
      # Number of points
  del.t <- 1/10
  sig.1 <- sig.2 <- .5
  q.1 <- q.2 <- 1
      # Model parameters
  FF <- diag(4)
  FF[1,3] <- FF[2,4] <- del.t
      # The location and velocity vectors
  GG <- matrix( 0, 2, 4 )
  GG[1,1] <- GG[2,2] <- 1
      # Unbiased measurement
  f <- function( x ) FF %*% x
  g <- function( x ) GG %*% x
      # The evolutions
  Q <- matrix( c( q.1 * del.t ^ 3 / 3, 0, q.1 * del.t ^ 2 / 2, 0,
                  0, q.2 * del.t ^ 3 / 3, 0, q.2 * del.t ^ 2 / 2,
                  q.1 * del.t ^ 2 / 2, 0, q.1 * del.t,         0,
                  0, q.2 * del.t ^ 2 / 2, 0,         q.2 * del.t ), 4, 4, byrow=T )   # State
  R <- matrix( c( sig.1 ^ 2, 0, 0, sig.2 ^ 2 ), 2 , 2)   # Obs
      # Variances
  x.0 <- c(10,-10,0,0)
  sim <- sim.data( x.0, f, g, R, Q, K )
      # Generate the data

  m.0 <- x.0
  P.0 <- Q
      # Initialize the mean and covariance
  base <- fkf( m.0, P.0, matrix(0,4,1), matrix(0,2,1), FF, GG, Q, R, sim$y )
      # From the FKF package
  alt <- ukf.compute( m.0, P.0, sim$y, f, g, Q, R, length(m.0), quad=F )
      # the unscented filter
  expect_equal( base$att, alt$m )
  expect_equal( base$at, alt$m.pred )
  expect_equal( base$Ptt, alt$P )
  expect_equal( base$Pt, alt$P.pred )
  
#   plot( t(sim$x[1:2,]), type='l' )
#   points( t(sim$y), pch=19, col='blue' )
#   lines( t(alt$m[1:2,]), col='red' )
#       ## Plot
})


test_that("A two-dimensional non-linear example",{
# The penduluum example from Ch5 of Sarkka
  
  K <- 500
  del.t <- .01
  grav <- 9.8
  q.c <- .01
  f <- function( x ) c( x[1] + del.t * x[2], x[2] - grav * sin(x[1]) * del.t )
  g <- function( x ) sin( x[1] )
  Q <- matrix( c( q.c * del.t ^ 3 / 3, q.c * del.t ^ 2 / 2, q.c * del.t ^ 2 / 2, q.c * del.t ), 2, 2 )
  R <- matrix( .1, 1, 1 )
  x.0 <- c( 2, 0 )
  sim <- sim.data( x.0, f, g, R, Q, K )
  m.0 <- .8 * x.0
  P.0 <- 200*Q
  ukf <- ukf.compute( m.0, P.0, sim$y, f, g, Q, R, length(m.0), alpha=1, quad = F )
  
  plot(del.t*1:K, sim$x[1,], type='l', lwd=2 )
  points(del.t*1:K, sim$y, pch=19, col='blue' )
  lines(del.t*1:K, ukf$m[1,], lty=2, col='red' )
  
  
  
})