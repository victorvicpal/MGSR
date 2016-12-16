coefcorr <- function(RES)
{
	p <- dim(RES$B[[1]])[1]
	fun <- length(RES$B)
	comb <- t(cbind(rbind(1:p,1:p),combn(1:p, 2)))
	ind <- comb[(p+1):(p*(p+1)/2),]
	Bvpri <- RES$Bvector[1:p,]
	Bvsec <- cbind(RES$Bvector[(p+1):(p*(p+1)/2),])
	v <- rep(0,p*(p+1)/2-p)

	coef <- function(index)
	{
	for (j in 1:length(Bvsec[,1])){
	v[j] <- Bvsec[j,index]/sqrt(abs(Bvpri[ind[j,1],index])*abs(Bvpri[ind[j,2],index]))
	}
	v
	}

	V <- lapply(1:fun,coef)

	comb <- cbind(comb,rep(0,length(comb[,1])))
	comb[1:p,3] <- 1

	reorder <- function(index)
	{
		CC <- list(0)
		s <- matrix(0,p,p)
		comb[(p+1):(p*(p+1)/2),3] <- V[[index]]
		s[comb[,1:2]] <- comb[,3]
		s[lower.tri(s)] = t(s)[lower.tri(s)]
		colnames(s) <- unique(unlist(strsplit(rownames(RES$Bvector)[1:p]," ")))
		rownames(s) <- colnames(s)
		CC <- s
		CC
	}

	CC <- lapply(1:fun,reorder)
	names(CC) <- names(RES$B)
	CC
}