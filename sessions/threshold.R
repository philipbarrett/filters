library(scales)

set.seed(4321)
sig.eps <- 1
rho <- .9
theta.hat <- 0

Q <- sig.eps^2
R <- .0 ^ 2
f <- function(x) rho * x
g <- function(x) if( x > theta.hat ) 1 else 0

x.0 <- 0
K <- 100

v.y <- v.x <- rep(0, K)
v.x[1] <- x.0
v.y[1] <- g( v.x[1] ) + rnorm( 1, 0, R )
for( i in 2:K ){
  v.x[i] <- f( v.x[i-1] ) + rnorm( 1, 0, sqrt( Q ) )
  v.y[i] <- g( v.x[i] ) + rnorm( 1, 0, sqrt( R ) )
}


kappa <- 10
thresh <- ukf.compute( mean(v.x), var(v.x), v.y, f, g, Q, R, 1, alpha=1, kappa=kappa, quad=F )

plot( c(1,K), range( c(thresh$m, v.x) ), type='n' )
points( 1:K, 1.5 * sd(v.x) * ( 2*v.y-1 ), pch=19, col=alpha('blue', .5), cex=.5 )
lines( 1:K, thresh$m.pred[-K], col='red', lwd=2 )
lines( 1:K, thresh$m + sqrt( c( thresh$P.pred[-K] ) ), col='red', lwd=2, lty=2 )
lines( 1:K, thresh$m - sqrt( c( thresh$P.pred[-K] ) ), col='red', lwd=2, lty=2 )
lines( 1:K, v.x, lwd=2 )
abline( h=0, lwd=1 )

# plot( 1:K, sqrt(c( thresh$P )), type='l' )
# points( 1:K, min(sqrt(c( thresh$P ))) + diff(range(sqrt(c( thresh$P )))) * v.y, pch=19, col='blue' )
