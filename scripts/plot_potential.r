par(mfrow=c(2,1))
min_x = -1
max_x = 1
min_y = -1
max_y = 1
x <- seq( min_x, max_x, len=100 )
y <- seq( min_y, max_y, len=100 )
z <- matrix( sin( 2*pi*(expand.grid(x,y)$Var1 + expand.grid(x,y)$Var2)) + cos( 2*pi*(expand.grid(x,y)$Var1 + expand.grid(x,y)$Var2)),length(x),length(y))
w <- exp(-z)
Lab.palette <- colorRampPalette(c("Black","Red", "Yellow"), space = "rgb")
cont<-filled.contour(x,y,w, xlab="x", ylab="y", nlevels=50, col=Lab.palette(128), main="Potential Exp(-V(x,y))")

library(hexbin)
endpoints = read.table("end_points.txt")
x <- ts(endpoints[1:1])
y <- ts(endpoints[2:2])

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))

h<-hexbinplot(x~y, xbins=500, xlim=c(-1,1), ylim=c(-1,1))


library(hexbin)

library(gridExtra)    
do.call(grid.arrange, c(c(h,cont), ncol=2))