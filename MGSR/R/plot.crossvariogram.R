plot.crossvariogram <- function(CV,RES=NULL)
{
	n <- length(CV)-3
	comb <- cbind(CV$comb,1:n)

	for (i in 1:length(comb[,1]))
	{
		if (comb[i,1]!=comb[i,2])
		{
			comb <- rbind(comb,c(rev(comb[i,1:2]),0))
		}
	}

	comb <- comb[order(comb[,1],comb[,2]),]

	myplot <- function(index)
	{
		if (index==0){frame()}
		else
		{
		plot(CV[[index]][1:2],main=names(CV)[index],col='blue')
		if(is.null(RES)==FALSE)	
			{
				lines(CV[[1]][,1],t(as.matrix(RES$Bvector[index,]))%*%t(as.matrix(RES$FUN)))
			}
		}
	}

	par(mfrow=c(CV$variables,CV$variables))
	par(mar=c(1.8,1.8,1.8,1.8))

	lapply(comb[,3],FUN=myplot)

}