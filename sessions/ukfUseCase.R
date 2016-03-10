
Q <- .25
R <- .25
f <- function(x) x - cos(x)
g <- function(x) x

x.0 <- -1.5
K <- 100

v.y <- v.x <- rep(0, K)
v.x[1] <- x.0
v.y[1] <- g( v.x[1] ) + rnorm( 1, 0, R )
for( i in 2:K ){
  v.x[i] <- f( v.x[i-1] ) + rnorm( 1, 0, Q )
  v.y[i] <- g( v.x[i] ) + rnorm( 1, 0, R )
}

plot( f, xlim=c(-2,2), ylim=c(0, max(v.x)) )
points( v.x[-K], v.x[-1] )
abline(0,1)

plot( 1:K, v.x, type='l' )
plot( v.x, v.y, pch=19, col='red' )

ukf <- ukf.compute( mean(v.x), var(v.x), v.y, f, g, Q, R, 1 )
plot( 1:K, v.x, type='l' )
points( 1:K, log( v.y ), pch=19, col='blue' )
lines( 1:K, ukf$m[,-1], col='red' )
