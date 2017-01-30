n_pairs_opt <- function(coord,values,n_min=12,n_max=20)
{
	v <- rep(0,(n_max-n_min+1))
	n <- 1
	for (i in n_min:n_max)
	{
		CV <- crossvariogram(coord,values,i)
		aux <- CV[[1]][,3]
		v[n] <- max(aux)-min(aux)
		n <- n+1
	}
	best_ap <- which.min(v)
	list(dif_pairs=v,best=best_ap+n_min-1)
}