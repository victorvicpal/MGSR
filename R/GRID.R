GRID <- function(DAT,lag=0.25,x1=0,x2=0,y1=0,y2=0)
{
	x.range <- (range(DAT[,1])) + c(x1,x2)
	y.range <- (range(DAT[,2])) + c(y1,y2)
	xygrid <- expand.grid(x=seq(from=x.range[1], to=x.range[2], by=lag),y=seq(from=y.range[1], to=y.range[2], by=lag))
}
