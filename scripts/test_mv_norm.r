library(MASS)

# The hexbin
sigma <- matrix( c(1,0,0,1),2,2)
samples <- mvrnorm(n = 5000, c(0,0), sigma)
x <- samples[,1]
y <- samples[,2]
hexbinplot(x~y, xbins=80, xlim=c(-3,3), ylim=c(-3,3))

# The background
x <- seq( -3, 3, len=100 )
y <- seq( -3, 3, len=100 )
w <- matrix( expand.grid(x,y)$Var1^2/2 + expand.grid(x,y)$Var2^2/2,length(x),length(y))
z <- 1/(sqrt(2*pi))*exp( -w )

filled.contour(x,y,z, xlab="x", ylab="y", nlevels=50, col=Lab.palette(128) )