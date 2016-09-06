par(mfrow=c(2,1))
min_x = -1
max_x = 1
min_y = -1
max_y = 1
x <- seq( min_x, max_x, len=100 )
y <- seq( min_y, max_y, len=100 )
z <- matrix( sin( 2*pi*(expand.grid(x,y)$Var1 + expand.grid(x,y)$Var2)) + cos( 2*pi*(expand.grid(x,y)$Var1 + expand.grid(x,y)$Var2)),length(x),length(y))

#Lab.palette <- colorRampPalette(c("Black","Red", "Yellow"), space = "rgb")
#cont<-filled.contour(x,y,z, xlab="x", ylab="y", nlevels=50, col=Lab.palette(128), main="V(x,y)")

w <- exp(-z)
Lab.palette <- colorRampPalette(c("Black","Red", "Yellow"), space = "rgb")
cont<-filled.contour(x,y,w, xlab="x", ylab="y", nlevels=90, col=Lab.palette(128), main="Exp(-V(x,y))")


library(hexbin)


plotEndPointsHex <- function( file_name )
{
  eps = read.table(file_name)
  x <- ts(eps[1:1])
  y <- ts(eps[2:2])
  hexbinplot(x~y, xbins=100, xlim=c(-1,1), ylim=c(-1,1))
}

list_of_files = list.files(pattern="end_points.*.txt")
lapply(list_of_files, plotEndPointsHex)
