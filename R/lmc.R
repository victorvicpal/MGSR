lmc <- function(CV,fun,a,tol=0.001,mode='aut')
{
	#Valores iniciales
	n <- length(CV$data[,1]) #Número de elementos medidos
	k <- length(CV[[1]][,1]) #Número de lags
	m <- length(fun)+1	#Número de funciones + Nugget
	p <- CV$variables	#Número de variables medidas
	pin <- (p*(1+p)/2)

	#Paquetes necesarios
	require(matrixcalc)
	#Ajustar b y nug de manera natural
	#RESULT <- GmodCross2(CV,fun,a,mode)
	RESULT <- GmodCross3(CV,fun,a)

	#Matrices A y Dist
	Dist <- as.matrix(dist(CV$data[,1:2]))
	A <- Amatrix(CV$data[,1:2],n=length(CV[[1]][,1]))

	#Calcular TRAZA 
	#Matriz kxk trazas inicial (depende de fun, a y A(h))
	##Combinación funciones y h
	CombFUN <- t(cbind(rbind(1:m,1:m),combn(1:m,2))) #Combinaciones entre funciones
	CombH <- t(cbind(rbind(1:k,1:k),combn(1:k,2)))	#Combinaciones entre lags
	CombFUNH <- expand.grid(1:m,1:k)	#Mezcla
	names(CombFUNH) <- c('FUN','H')
	##Matrices nxn de las funciones "tipo" 
	###OJO! Si hay más de una función repetida hay que tener en cuenta el rango de la misma
	#Cambios en las funciones el 20/05/2016 En el github sigue la versión antigua.
	CorrSph <- function(D,as){Fsph <- (ifelse(D>as,Fsph <- 0,Fsph <- (1-1.5*(D/as)+0.5*(D/as)^3)))}#MIRAR SI ES 1-1.5...
	CorrExp <- function(D,as){Fexp <- (exp(-D/as))}
	CorrGau <- function(D,as){Fgau <- (exp(-(D/as)^2))}
	CorrLin <- function(D,as){Flin <- (ifelse(D>as,Flin <- 0,Flin <- (1-(D/as))))} #lineal
	CorrPow <- function(D,as){Fpow <- (ifelse(as<=2 & as>0,Fsph <- (1-D^as),stop('el rango no se encuentra entre 0 y 2')))}
	##FUN
	FUN <- list(diag(1,n,n))
	for (i in 1:length(fun))
	{if(fun[i]=='Exp'){FUN[[i+1]]<-CorrExp(Dist,a[i])}
	else{if(fun[i]=='Sph'){FUN[[i+1]]<-CorrSph(Dist,a[i])}else{FUN[[i+1]]<-CorrGau(Dist,a[i])}}}
	names(FUN) <- c('Nug',fun)

	prevtraza <- function(v1,v2)
	{list(FUN[[v1]]%*%A[[v2]])}
	
	FUNH <- mapply(prevtraza,CombFUNH[,1],CombFUNH[,2])

	traza <- function(w1)
	{
		ind1 <- rep(which(CombFUNH$H==CombH[w1,1]), each=length(fun))
		ind2 <- rep(which(CombFUNH$H==CombH[w1,2]),each=length(fun))
		r <- matrix.trace(FUNH[[ind1[1]]]%*%FUNH[[ind2[1]]])
		for (i in 2:length(ind1))
		{r <- r+matrix.trace(FUNH[[ind1[i]]]%*%FUNH[[ind2[i]]])}
		t <- r
	}

	CombH <- cbind(CombH,unlist(lapply(1:length(CombH[,1]),traza)))

	TRAZAS <- matrix(0,k,k)
	TRAZAS[lower.tri(TRAZAS,diag=TRUE)] <- CombH[,3]
	TRAZAS[upper.tri(TRAZAS)] = t(TRAZAS)[upper.tri(TRAZAS)] #Matriz kxk TRAZAS

	#Calcular la primeras matrices de B0 a partir de RESULT
	comb <- t(cbind(rbind(1:p,1:p),combn(1:p, 2)))
	comb <- cbind(comb,rep(0,(p*(1+p)/2)))

	Bnew <- RESULT$B
	Bnewlast <- Bnew

	bcalc <- function(B)
	{sum((Reduce('+',B))^2)} #Cálculo de la constante b que multiplica a la matriz de trazas

	b <- bcalc(Bnew)

	WSS <- function(index)
	{bg <- Bnew[index,1]*RESULT$FUN[,1]
		for (i in 2:m)
		{bg <- bg+Bnew[index,i]*RESULT$FUN[,i]} 
	wss <- t(CV[[index]][,2]-bg)%*%solve(b*TRAZAS)%*%(CV[[index]][,2]-bg)#ANTES TENIA bg[[index]]
	}

	SUMWSS <- function(comb,wss)
	{
		s <- matrix(0,p,p)
		comb[,3] <- wss
		s[comb[,1:2]] <- comb[,3]
		s[lower.tri(s)] = t(s)[lower.tri(s)]
		sum(s)
	}

	wss0 <- unlist(lapply(1:pin,WSS))
	WSS1 <- SUMWSS(comb,wss0)
	WSSInic <- WSS1
	WSS0 <- 0

	getbg <- function(index)
	{	
		bg <- Bnewlast[index,1]*RESULT$FUN[,1]
		for (i in 2:m)
		{bg <- rbind(bg,Bnewlast[index,i]*RESULT$FUN[,i])}
		rownames(bg) <- c('Nug',fun)
		t(bg)
	}

	getgammas0 <- function(index, funi)
		{bg <- lapply(1:pin,getbg)
		bgd <- CV[[index]][,2]-(rowSums(bg[[index]])-bg[[index]][,funi])
		bgd}

	fitmodel <- function(index,funi)
	{
	gammas0 <- mapply(getgammas0,1:pin,funi)
	bnew <- solve(t(RESULT$FUN[,funi])%*%solve(b*TRAZAS)%*%RESULT$FUN[,funi])%*%t(RESULT$FUN[,funi])%*%solve(b*TRAZAS)%*%gammas0[,index]
	}

	Bt <- function(comb,bfit)
	{
		comb[,3] <- bfit
		s <- matrix(0,p,p)
		s[comb[,1:2]] <- comb[,3]
		s[lower.tri(s)] = t(s)[lower.tri(s)]
		s
	}

	DVS <- function(BN)
	{PC <- eigen(BN)
	PC$values[PC$values<0] <- 0
	ord <- order(PC$values,decreasing = T)
	PC$values <- PC$values[ord]
	PC$vectors <- PC$vectors[,ord]
	dG <- PC$vectors%*%diag(PC$values)%*%t(PC$vectors)}
	#Inversa generalizada

	Bvector <- function(B)
	{bp <- c(diag(B),B[lower.tri(B)])
	bp}

	it <- 0

	while (abs(WSS1-WSS0)>tol)
	{
		s <- 1
		while (s<=m)
		{
			bfit <- mapply(fitmodel,1:pin,s)
			Bfit <- Bt(comb,bfit)
			Bdvs <- DVS(Bfit)
			Bnew[,s] <- Bvector(Bdvs)
			b <- bcalc(Bnew)
			s <- s+1
		}
		wss1 <- unlist(lapply(1:(p*(1+p)/2),WSS))
		WSS2 <- SUMWSS(comb,wss1)
		if (abs(WSS1-WSS0) < abs(WSS2-WSS1)) break
		if (abs(WSS2)>abs(WSS1)) break
		WSS0 <- WSS1 #EL FALLO ESTÁ EN EL CONTROL DE WSS OJO!
		WSS1 <- WSS2
		Bnewlast <- Bnew
		it <- it+1
		#El problema es que WSS1-WSS0 puede ser mayor que el anterior, entonces no para y se aleja
	}
	
	orderBfinal <- function(comb,BT)
	{
		B <- list(0)
		s <- matrix(0,p,p)
	for (i in 1:m)
		{
			comb[,3] <- BT[,i]
			s[comb[,1:2]] <- comb[,3]
			s[lower.tri(s)] = t(s)[lower.tri(s)]
			B[[i]] <- s
		}
		names(B) <- c('Nug',fun)
		B
	}
	Bfinal <- orderBfinal(comb,Bnewlast)

	if(it==0)
	{
		for (i in 1:m)
		{
			Bfinal[[i]] <- DVS(Bfinal[[i]])
		}
		for (i in 1:m)
		{
			Bnewlast[,i] <- c(diag(Bfinal[[i]]),Bfinal[[i]][lower.tri(Bfinal[[i]])])
		}
	}

	NEWRESULT <- list(var=CV$variables,rango=a,B=Bfinal,Bvector=Bnewlast,FUN=RESULT$FUN,WSSini=WSSInic,WSSfinal=WSS1,iter=it,data=CV$data)

	#if( NEWRESULT$WSSini<NEWRESULT$WSSfinal ) stop('Valores iniciales mal elegidos. Prueba con nuevas combinaciones')

	NEWRESULT
}